docker build -t orph-env .

docker run --rm \
  -p 8787:8787 \
  -e PASSWORD=mypass \
  -v "$(pwd)":/home/rstudio/work \
  orph-env

  #After running the above command, you can access RStudio Server by navigating to http://localhost:8787 in your web browser. Use "rstudio" as the username and "mypass" as the password to log in.