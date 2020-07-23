# Reverting development files after a test gone wrong

If reverting the files in a data delivery test fails,
then the files delivered by the test
will no longer be in the broad-gotc-dev-storage bucket.
The `revert.sh` script can be used
to copy prod data to the dev bucket
to restore the files.
The firecloud bucket and workspace created by the test
will still need to be deleted manually.
