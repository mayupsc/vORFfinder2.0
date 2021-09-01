docker stop shiny 
docker rm shiny
docker run --name shiny -d -p 8383:8383 --restart=always -v ${PWD}:/srv/shiny-server/ yumapsc/vorffinder2.1
#sudo docker exec -ti shiny /bin/bash
