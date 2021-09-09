# Set git config information
git config --global user.name "James Balamuta"
git config --global user.email "balamut2@illinois.edu"

# Clone the gh-pages repository
git clone -b gh-pages \
  https://${GH_TOKEN}@github.com/smac-group/docs.git \
  doc-output
  
# Change to the gh-page clone doc-output directory
cd doc-output/gmwm

# Copy generated output to doc-output
cp -r ../../docs/* ./

# Add all files to the repo
git add *
git commit -a -m "Updating gmwm-docs (${TRAVIS_BUILD_NUMBER})"
git push origin gh-pages