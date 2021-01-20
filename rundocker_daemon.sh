docker stop shiny2 
docker rm shiny2 
docker run --name shiny2 -d -p 8383:8383 -v $PWD:/srv/shiny-server/ yumapsc/vorffinder2.1
sudo docker exec -ti shiny2 /bin/bash
