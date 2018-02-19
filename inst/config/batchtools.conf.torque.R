cluster.functions = makeClusterFunctionsTORQUE(
                      template = pulsar::findTemplateFile('simpletorque'),
                      scheduler.latency = 1,
                      fs.latency = 65)
