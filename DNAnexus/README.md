### Guides for creating DNAnexus apps.


Links:



dx-toolkit commands:

Get input for app using [jq](https://stedolan.github.io/jq/) library for bash
```bash
dx describe pbgermline --json | jq -r .inputSpec > parabricks.germline.inputTemplate.json
```
