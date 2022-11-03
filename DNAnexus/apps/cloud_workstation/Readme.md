# Cloud Workstation App

## What does this app do?
This app sets up a cloud workstation which you can access by running the app with the **--ssh** or **--allow-ssh** flags.

## What are typical use cases for this app?
This app can be used as a workstation inside of the DNAnexus cloud platform.  By running the app with **--ssh** or
**--allow-ssh**, users can login to a machine inside of the DNAnexus cloud platform.  From there, users can upload/download
data to/from the project in which the app is run, perform data analysis, and install additional packages from sources such as 
apt, cran, pip, github, etc.

Note: in order to access files stored in the project in which this app is being run, users must provide the project-identifier.
This can easily be done with a special environment variable called $DX\_PROJECT\_CONTEXT\_ID.  For instance, to download a file called
myreads.fastq.gz from the parent project, users would simply run 
```dx download $DX_PROJECT_CONTEXT_ID:/myreads.fastq.gz```

## What are the inputs?
The user can provide a maximum session length value.  After this amount of time has passed, the workstation will automatically 
shut-down.  Timeout is provided using suffixes s, m, h, d, w, M, y.

During a session, users can check how much time remains until the session times out by running
```dx-get-timeout```
Users can reset the timeout by running
```dx-set-timeout <timeout>```
using the suffixes s, m, h, d, w, M, y.

Additionally, users can provide a list of files to download at the app startup.  These files will be downloaded to the home
directory automatically, allowing easy access for data analysis.

Finally, users can create snapshots of their cloud environment by running
```dx-create-snapshot <snapshot name (optional)>```
This call will generate snapshot files, which can be provided as input to future runs of cloud_workstation to pick up where you left off.

Note: `dx-snapshot` has been deprecated. Please use `dx-create-snapshot` instead.
## Examples

### Launch cloud_workstation with SSH enabled, set to run for 10 hours
```
$ dx run cloud_workstation -imax_session_length=10h --allow-ssh --brief -y --name "Alice's workstation"
job-xxxx
```

After the worker is launched and is running, SSH into it:
```
$ dx ssh job-xxxx
```

Alternatively, you can use `--ssh` to directly SSH into the worker when launching it. So the previous two steps are merged into:
```
$ dx run cloud_workstation -imax_session_length=10h --ssh --brief -y --name "Alice's workstation"
```


### When SSH-ed in the cloud_workstation worker

Check time left for the worker:
```
$ dx-get-timeout
0 days 9 hours 27 minutes 44 seconds
```

To update the time period for the worker to be up:
```
$ dx-set-timeout 20h 
$ dx-get-timeout
0 days 19 hours 59 minutes 55 seconds
```

To save the entire workspace on the worker as a cloud_workstation snapshot:
```
$ dx-create-snapshot 
/dx_stderr
/dx_stdout
/etc
/etc/mtab
/etc/resolv.conf
[... more files ...]
```

You can then find the snapshot file saved at the root level of the project in which the job is run and named with a timestamp (e.g. Dec_20_2021_19_01.snapshot).
You can also save the snapshot with a custom file name using option `--name` (`dx-create-snapshot --name "snapshot_2021"`) and it will be saved as "snapshot_2021".

To terminate and close the worker before the set timeout:
```
$ dx-set-timeout 0s
```

### Re-use a cloud_workstation snapshot

To continue working on a previously created cloud_workstation snapshot:
```
$ dx run cloud_workstation -isnapshot="Dec_20_2021_19_01.snapshot" -imax_session_length=10h --allow-ssh --brief -y --name "Alice's workstation"
```
This will launch a cloud_workstation for you to SSH in, download the snapshot from previous run, and decompress the workspace to where it was.
The details field of the snapshot file description should contain the semantic version (i.e. 2.0.0) of the cloud_workstation that created the snapshot. A snapshot may not be compatible with a different major version of the cloud_workstation, but you can use the snapshot with the specific version X.Y.Z of the cloud_workstation that created it by running `dx run app-cloud_workstation/X.Y.Z`.

