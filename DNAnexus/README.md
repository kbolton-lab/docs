### Guides for creating DNAnexus apps.


Links:



**dx-toolkit commands:**

Get input for app using [jq](https://stedolan.github.io/jq/) library for bash
```bash
dx describe pbgermline --json | jq -r .inputSpec > parabricks.germline.inputTemplate.json
```

a. dx-docker -h
```bash
usage: dx-docker [-h] [-q] {pull,run,add-to-applet,create-asset} ...

positional arguments:
  {pull,run,add-to-applet,create-asset}
    pull                Pulls a docker image for use in DNAnexus
    run                 Runs a docker image in a container
    add-to-applet       Adds a local Docker image to an applet
    create-asset        Caches a local Docker image as an asset in the DNAnexus platform (subject to change)

optional arguments:
  -h, --help            show this help message and exit
  -q, --quiet           Suppress printing pull progress to stderr
  ```
