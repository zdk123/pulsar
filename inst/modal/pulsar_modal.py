"""Modal future backend helper for pulsar (see ?pulsar::modal).

R starts this file as a long-lived process:

    python pulsar_modal.py serve

It opens an ephemeral Modal app for the R session and reads requests from
stdin, one per line, tab-separated:

    submit <key> <inpath> <outpath>   ->  ok <call id>   |  error <message>
    cancel <key>                      ->  ok
    quit

`inpath` holds a serialized future (its data from future::getFutureData()).
Each submitted future runs `pulsar:::.modal_eval()` in its own container.
When it finishes, the helper writes the FutureResult to `outpath`, or the
failure message to `outpath` + ".err" when the container itself failed (lost
worker, timeout, out of memory). Files are written under a temporary name and
renamed, so R sees each one complete or not at all.

The worker image and resources are set by PULSAR_MODAL_* environment
variables, which the R side fills from plan(pulsar::modal, ...).
"""

import contextlib
import os
import pathlib
import random
import re
import subprocess
import sys
import tempfile
import threading
import time

import modal

ENV = os.environ.get


def _lines(name, default=""):
    return [x for x in ENV(name, default).split("\n") if x.strip()]


R_VERSION = ENV("PULSAR_MODAL_R_VERSION", "4.6.2")
# package specs, one per line: "pkg" (CRAN), "bioc::pkg" (Bioconductor) or
# "owner/repo@ref" (GitHub)
R_PACKAGES = _lines("PULSAR_MODAL_R_PACKAGES", "huge")
# local package directories, one per line, in install (dependency) order
LOCAL = _lines("PULSAR_MODAL_LOCAL")
# system libraries; libglpk is needed by igraph, which huge imports
APT = [p for p in ENV("PULSAR_MODAL_APT", "libglpk-dev").split(",") if p]
# Modal secrets available to the package installs (e.g. one holding GITHUB_PAT)
SECRETS = [s for s in ENV("PULSAR_MODAL_SECRETS", "").split(",") if s]
BASE_IMAGE = ENV("PULSAR_MODAL_IMAGE", "")  # registry image with R, replaces rocker
DOCKERFILE = ENV("PULSAR_MODAL_DOCKERFILE", "")  # or a Dockerfile building one
PULSAR_SOURCE = ENV("PULSAR_MODAL_SOURCE", "")  # local package directory
PULSAR_REF = ENV("PULSAR_MODAL_REF", "")  # GitHub ref, e.g. a commit SHA
CPU = float(ENV("PULSAR_MODAL_CPU", "1"))
MEMORY = int(ENV("PULSAR_MODAL_MEMORY", "2048"))
TIMEOUT = int(ENV("PULSAR_MODAL_TIMEOUT", "3600"))
RETRIES = int(ENV("PULSAR_MODAL_RETRIES", "0"))
MAX_CONTAINERS = int(ENV("PULSAR_MODAL_MAX_CONTAINERS", "100"))
APP_NAME = ENV("PULSAR_MODAL_APP", "pulsar")

# Specs are pasted into R code below, so only these shapes are accepted (the
# R side checks them too, with the same patterns)
_NAME = r"[A-Za-z][A-Za-z0-9.]*"
SPEC_CRAN = re.compile(rf"^{_NAME}$")
SPEC_BIOC = re.compile(rf"^bioc::({_NAME})$")
SPEC_GITHUB = re.compile(r"^[A-Za-z0-9][A-Za-z0-9-]*/[A-Za-z0-9_-][A-Za-z0-9_.-]*"
                         r"(/[A-Za-z0-9_][A-Za-z0-9_./-]*)?(@[A-Za-z0-9_][A-Za-z0-9_./-]*)?$")
SPEC_REF = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_./-]*$")


def _rvec(xs):
    return "c(" + ",".join(f"'{x}'" for x in xs) + ")"


def _split_specs(specs):
    cran, bioc, github = [], [], []
    for x in specs:
        if ".." in x:
            raise ValueError(f"package spec may not contain '..': {x!r}")
        if SPEC_CRAN.match(x):
            cran.append(x)
        elif SPEC_BIOC.match(x):
            bioc.append(SPEC_BIOC.match(x).group(1))
        elif SPEC_GITHUB.match(x):
            github.append(x)
        else:
            raise ValueError(f"not a CRAN, bioc:: or owner/repo@ref package spec: {x!r}")
    return cran, bioc, github


def _rscript(code):
    return f"Rscript -e \"{code}\""


cran, bioc, github = _split_specs(R_PACKAGES)
if PULSAR_REF and (not SPEC_REF.match(PULSAR_REF) or ".." in PULSAR_REF):
    raise ValueError(f"bad pulsar ref {PULSAR_REF!r}")
secrets = [modal.Secret.from_name(s) for s in SECRETS]

if DOCKERFILE:
    image = modal.Image.from_dockerfile(DOCKERFILE, add_python="3.12")
else:
    image = modal.Image.from_registry(BASE_IMAGE or f"rocker/r-ver:{R_VERSION}",
                                      add_python="3.12")
if APT:
    image = image.apt_install(*APT)
# BiocManager always: GitHub and local packages may import Bioconductor ones
deps = ["remotes", "Matrix", "future", "BiocManager"] + cran
# repos for remotes installs: CRAN (the image's Package Manager) and Bioconductor
_REPOS = "options(repos=BiocManager::repositories()); "
image = image.run_commands(
    _rscript(f"install.packages({_rvec(deps)})"),
    _rscript(f"stopifnot(all(sapply({_rvec(deps)}, requireNamespace)))"),
)
if bioc:
    image = image.run_commands(
        _rscript(f"BiocManager::install({_rvec(bioc)}, update=FALSE, ask=FALSE); "
                 f"stopifnot(all(sapply({_rvec(bioc)}, requireNamespace)))"),
        secrets=secrets,
    )

# not uploaded from package directories (R/modal.R's .modalfiles skips the
# same files): VCS and check output, editor and renv state, and objects
# compiled on the controller, which would not link in the Linux image
_IGNORE = ["*.Rcheck", "*.tar.gz", ".git", "doc", "Meta", ".Rproj.user", "renv",
           ".Rhistory", "**/.DS_Store", "**/__pycache__", "**/*.o", "**/*.so", "**/*.dll"]
if PULSAR_SOURCE:
    image = image.add_local_dir(PULSAR_SOURCE, "/opt/pulsar-src", copy=True,
                                ignore=_IGNORE
    ).run_commands("R CMD INSTALL --no-docs --no-multiarch /opt/pulsar-src")
elif PULSAR_REF:
    image = image.run_commands(
        _rscript(f"remotes::install_github('zdk123/pulsar@{PULSAR_REF}', "
                 "dependencies=FALSE, upgrade='never')"),
        secrets=secrets,
    )
elif modal.is_local():
    # Modal re-imports this file in each container without these variables;
    # the image is already built by then, so only the controller checks them
    raise RuntimeError("set PULSAR_MODAL_SOURCE or PULSAR_MODAL_REF")

# GitHub and local packages come after pulsar, so they may depend on it; local
# ones come last, so they may depend on GitHub ones (not the reverse). Their
# CRAN and Bioconductor dependencies (Depends, Imports, LinkingTo) are
# installed with them. A branch ref is cached at its
# first build; use a commit SHA, or rebuild=TRUE.
for spec in github:
    image = image.run_commands(
        # remotes stops (failing the build) if the install fails
        _rscript(_REPOS + f"remotes::install_github('{spec}', upgrade='never', dependencies=NA)"),
        secrets=secrets,
    )
# a local package's files are part of the image, so changing one rebuilds
# from there on
for i, path in enumerate(LOCAL):
    image = image.add_local_dir(path, f"/opt/local-src/{i}", copy=True,
                                ignore=_IGNORE
    ).run_commands(
        _rscript(_REPOS + f"remotes::install_local('/opt/local-src/{i}', upgrade='never', "
                 "build=FALSE, dependencies=NA)"),
        secrets=secrets,
    )

app = modal.App(APP_NAME, image=image)


@app.function(
    cpu=CPU,
    memory=MEMORY,
    timeout=TIMEOUT,
    retries=RETRIES,
    max_containers=MAX_CONTAINERS,
)
def run_future(payload: bytes) -> bytes:
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        (d / "in.rds").write_bytes(payload)
        proc = subprocess.run(
            ["Rscript", "-e", "a <- commandArgs(TRUE); pulsar:::.modal_eval(a[1], a[2])",
             str(d / "in.rds"), str(d / "out.rds")],
            capture_output=True,
            text=True,
        )
        out = d / "out.rds"
        if proc.returncode != 0 or not out.exists():
            raise RuntimeError(
                f"R worker exited with status {proc.returncode}: "
                + (proc.stderr.strip()[-2000:] or "no output")
            )
        return out.read_bytes()


def _write_atomic(path: pathlib.Path, data: bytes):
    tmp = path.with_name(path.name + ".part")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _collect(call, outpath: pathlib.Path, done):
    try:
        _write_atomic(outpath, call.get())
    except BaseException as e:  # noqa: BLE001 - every failure goes back to R
        msg = f"{type(e).__name__}: {e}".strip()
        _write_atomic(outpath.with_name(outpath.name + ".err"),
                      (msg or "modal job failed").encode())
    finally:
        done()


def _start_app(stack, attempts=8):
    # Modal limits how fast an account creates apps, so sessions started
    # together (one per simulation draw, say) can be refused; wait and retry
    for i in range(attempts):
        try:
            return stack.enter_context(app.run())
        except modal.exception.ResourceExhaustedError as e:
            if i == attempts - 1:
                raise
            wait = min(60, 2 ** i) * random.uniform(0.5, 1.5)
            print(f"pulsar: {e} (retrying in {wait:.0f}s)", file=sys.stderr, flush=True)
            time.sleep(wait)


def serve():
    # stdout carries the protocol; send everything else (Modal's progress
    # output, prints) to stderr
    proto = os.fdopen(os.dup(1), "w", buffering=1)
    os.dup2(2, 1)
    sys.stdout = sys.stderr

    def reply(*fields):
        proto.write("\t".join(str(f).replace("\t", " ").replace("\n", " ")
                              for f in fields) + "\n")
        proto.flush()

    calls = {}
    with contextlib.ExitStack() as stack:
        stack.enter_context(modal.enable_output())
        _start_app(stack)
        reply("ready")
        for line in sys.stdin:
            req = line.rstrip("\n").split("\t")
            op = req[0]
            try:
                if op == "submit":
                    key, inpath, outpath = req[1], req[2], pathlib.Path(req[3])
                    call = run_future.spawn(pathlib.Path(inpath).read_bytes())
                    calls[key] = call
                    threading.Thread(target=_collect, daemon=True,
                                     args=(call, outpath,
                                           lambda k=key: calls.pop(k, None))).start()
                    reply("ok", call.object_id)
                elif op == "cancel":
                    call = calls.pop(req[1], None)
                    if call is not None:
                        call.cancel(terminate_containers=True)
                    reply("ok")
                elif op == "quit":
                    reply("ok")
                    break
                else:
                    reply("error", f"unknown request {op!r}")
            except Exception as e:  # noqa: BLE001
                reply("error", f"{type(e).__name__}: {e}")


if __name__ == "__main__" and sys.argv[1:] == ["serve"]:
    serve()
