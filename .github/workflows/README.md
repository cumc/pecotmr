## How to tag new releases with the "Upload new release" GitHub Action
When we are ready to tag a new release, we use this action to create in the repository and upload the source archive as release available on the [Releases](https://github.com/cumc/pecotmr/releases) page.
1. Make a pull request to update the version in the DESCRIPTION file.  Do **not** include any dashes in the version number! Do **not** manually tag a release for this version as it will be done later automatically by the workflow itself. 
2. Go to "Actions" at the top of the repository page.
3. Go to the "Upload new release" action on the side bar.
4. Click on the "Run workflow" drop down menu in the blue highlighted area
5. In the drop down menu, fill in the version number that you want to tag.  Make sure it is the same version as that in the DESCRIPTION file.
   - The commit checksum is optional --- by default it will use the latest commit to repository, but you can specify an older commit if needed.  Do **not** use a commit that is older than the commit used for the current version.
6. Click the green "Run workflow" button to dispatch the workflow.
## How to build a new conda package with the "Build conda package" GitHub Action
When we have tagged a new release, we use this action to build a new conda package.  The conda packages are currently uploaded to the [personal channel](https://anaconda.org/dnachun) of Daniel Nachun.  Upon submission to Bioconductor, a recipe will be submitted to bioconda and this workflow will be replaced by one to submit releases to Bioconductor, and bioconda will automatically update the package.
1. Use the "Upload new release" GitHub Action to tag a new release.  This workflow will fail if you try to build a package for a version which is not already tagged.
2. Go to "Actions" at the top of the repository page.
3. Go to the "Build conda package" action on the side bar.
4. Click on the "Run workflow" drop down menu in the blue highlighted area
5. In the drop down menu, fill in the version number that you want to build a package for.
   - The build version is optional and defaults to 0.  This can be incremented if the conda recipe has been changed but the tagged release has not been changed.  Do **not** change this setting without asking first - in most cases a new version needs to be tagged!
6. Click the green "Run workflow" button to dispatch the workflow.
