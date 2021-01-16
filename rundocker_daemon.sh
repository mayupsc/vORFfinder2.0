docker stop shiny2 
docker rm shiny2 
docker run --name shiny2 -d -p 3838:3838 -v /Users/mayu/Documents/huangtan/vORFfinder2.0/:/srv/shiny-server/ yumapsc/vorffinder2.1
sudo docker exec -ti shiny2 /bin/bash
