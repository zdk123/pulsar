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

import os
import pathlib
import subprocess
import sys
import tempfile
import threading

import modal

ENV = os.environ.get
R_VERSION = ENV("PULSAR_MODAL_R_VERSION", "4.6.2")
R_PACKAGES = [p for p in ENV("PULSAR_MODAL_R_PACKAGES", "huge").split(",") if p]
# system libraries; libglpk is needed by igraph, which huge imports
APT = [p for p in ENV("PULSAR_MODAL_APT", "libglpk-dev").split(",") if p]
PULSAR_SOURCE = ENV("PULSAR_MODAL_SOURCE", "")  # local package directory
PULSAR_REF = ENV("PULSAR_MODAL_REF", "")  # GitHub ref, e.g. a commit SHA
CPU = float(ENV("PULSAR_MODAL_CPU", "1"))
MEMORY = int(ENV("PULSAR_MODAL_MEMORY", "2048"))
TIMEOUT = int(ENV("PULSAR_MODAL_TIMEOUT", "3600"))
RETRIES = int(ENV("PULSAR_MODAL_RETRIES", "0"))
MAX_CONTAINERS = int(ENV("PULSAR_MODAL_MAX_CONTAINERS", "100"))
APP_NAME = ENV("PULSAR_MODAL_APP", "pulsar")


def _rvec(xs):
    return "c(" + ",".join(f"'{x}'" for x in xs) + ")"


deps = ["remotes", "Matrix", "future"] + R_PACKAGES
image = modal.Image.from_registry(f"rocker/r-ver:{R_VERSION}", add_python="3.12")
if APT:
    image = image.apt_install(*APT)
image = image.run_commands(
    f"Rscript -e \"install.packages({_rvec(deps)})\"",
    f"Rscript -e \"stopifnot(all(sapply({_rvec(deps)}, requireNamespace)))\"",
)
if PULSAR_SOURCE:
    image = image.add_local_dir(
        PULSAR_SOURCE,
        "/opt/pulsar-src",
        copy=True,
        ignore=["*.Rcheck", "*.tar.gz", ".git", "doc", "Meta", "**/__pycache__"],
    ).run_commands("R CMD INSTALL --no-docs --no-multiarch /opt/pulsar-src")
elif PULSAR_REF:
    image = image.run_commands(
        f"Rscript -e \"remotes::install_github('zdk123/pulsar@{PULSAR_REF}', "
        "dependencies=FALSE, upgrade='never')\""
    )
elif modal.is_local():
    # Modal re-imports this file in each container without these variables;
    # the image is already built by then, so only the controller checks them
    raise RuntimeError("set PULSAR_MODAL_SOURCE or PULSAR_MODAL_REF")

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
    with modal.enable_output(), app.run():
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
