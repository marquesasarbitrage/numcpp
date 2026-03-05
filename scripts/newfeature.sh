git checkout main && git pull origin main
git checkout -b "feature/$1"
git push -u origin "feature/$1"
echo "Branch feature/$1 created and pushed."