# cloud_workstation Developer Readme

<!--
TODO: Please edit this Readme.developer.md file to include information
for developers or advanced users, for example:

* Information about app internals and implementation details
* How to report bugs or contribute to development
-->

## Running cloud_workstation app with additional computational resources

This app has the following entry points:

* main

When running this app, you can override the instance type to be used by
providing the ``systemRequirements`` field to ```/applet-XXXX/run``` or
```/app-XXXX/run```, as follows:

    {
      systemRequirements: {
        "main": {"instanceType": "mem2_hdd2_x2"}
      },
      [...]
    }

See the [Run Specification](https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications#run-specification)
page in the API documentation for more information about the available instance types.

## Updating cloud_workstation dxapp.json and source code

If you are planning to create a new version of the app, please conform to [App(let) style guide](https://documentation.dnanexus.com/developer/apps/third-party-app-style-guide). If the changes might make older snapshots incompatible with the new version of cloud_workstation, for example, when the operating system on which the app is run is upgraded, you must increase the [MAJOR](https://semver.org/) version in dxapp.json (e.g., from 1.1.1 to 2.0.0).