## How to tag new releases with the "Upload new release" GitHub Action
When we are ready to tag a new release, we use this action to create in the repository and upload the source archive as release available on the [Releases](https://github.com/StatFunGen/pecotmr/releases) page.
1. Go to "Actions" at the top of the repository page.
2. Go to the "Upload new release" action on the side bar.
3. Click on the "Run workflow" drop down menu in the blue highlighted area
4. In the drop down menu, you do not need to modify any values if you just want to increment the patch version automatically.
    - We use semantic versioning of the form version X.Y.Z, where X is the major version, Y is the minor version and Z is the patch version
    - By default, the patch version is incremented because `increase_patch_version` is set to true
    - If you need to increment the minor or major version, this can be done by setting `increase_minor_version` or `increase_major_version` to true.  Please ask before changing these values.
    - A custom version can be specified instead, but please ask before using this.
    - In all cases the DESCRIPTION file will be updated automatically to the right version - do not manually change this file without asking.
    - The commit checksum is optional --- by default it will use the latest commit to repository, but you can specify an older commit if needed.  Do **not** use a commit that is older than the commit used for the current version.
5. Click the green "Run workflow" button to dispatch the workflow.
## How to build a new conda package with the "Build conda package" GitHub Action
When we have tagged a new release, we use this action to build a new conda package.  The conda packages are currently uploaded to the [personal channel](https://anaconda.org/dnachun) of Daniel Nachun.  Upon submission to Bioconductor, a recipe will be submitted to bioconda and this workflow will be replaced by one to submit releases to Bioconductor, and bioconda will automatically update the package.
1. Use the "Upload new release" GitHub Action to tag a new release.  This workflow will fail if you try to build a package for a version which is not already tagged.
2. Go to "Actions" at the top of the repository page.
3. Go to the "Build conda package" action on the side bar.
4. Click on the "Run workflow" drop down menu in the blue highlighted area
5. In the drop down menu, you do not need to modify any values if you just want to build the latest version of package.  If you need to build an older version of the package, you can enter a custom version, but please ask before doing this.
    - The build version is optional and defaults to 0.  This can be incremented if the conda recipe has been changed but the tagged release has not been changed.  Do **not** change this setting without asking first - in most cases a new version needs to be tagged!
6. Click the green "Run workflow" button to dispatch the workflow.
