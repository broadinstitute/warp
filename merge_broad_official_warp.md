# Step 1: From your project repository, check out a new branch and test the changes.

git checkout -b broadinstitute-develop develop
git pull https://github.com/broadinstitute/warp.git develop

# Step 2: Merge the changes and update on GitHub.

git checkout develop
git merge --no-ff broadinstitute-develop
git push origin develop