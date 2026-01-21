docker build -t orph-env .

docker run --rm \
  -p 8787:8787 \
  -e PASSWORD=mypass \
  -v $(pwd):/home/rstudio/work \
  orph-env