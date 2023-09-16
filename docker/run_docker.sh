#!/bin/bash
# Get the directory where the script resides
script_dir=$(dirname "${BASH_SOURCE[0]}")

# Change the working directory to the script's directory
cd "$script_dir/.."
echo $PWD

SF_NAME="test_imu"

# Define the name of the ROS Noetic container
container_name="noetic_srrg"
image_name="noetic_srrg"

if [ "$1" = "reset" ]; then
    # Check if a container with the specified name is running
    if [ "$(docker ps -q -f name=$container_name)" ]; then
        # A container with the specified name is running; stop and remove it
        echo "Stopping and removing the existing container $container_name..."
        docker stop $container_name
        docker rm $container_name
    elif [ "$(docker ps -aq -f status=exited -f name=$container_name)" ]; then
        # A container with the specified name exists but is stopped; remove it
        echo "Removing the existing stopped container $container_name..."
        docker rm $container_name
    fi
fi
    # Check if a container with the specified name is already running
    if [ "$(docker ps -q -f name=$container_name)" ]; then
        # A container with the specified name is running
        echo "Container $container_name is already running."
        echo "Attaching to its shell..."
        docker exec -it $container_name /bin/bash
    else
        # No container with the specified name is running
        # Check if a container with the specified name exists (stopped)
        if [ "$(docker ps -aq -f status=exited -f name=$container_name)" ]; then
            # A container with the specified name exists but is stopped
            echo "Starting existing container $container_name..."
            docker start $container_name 
            docker exec -it $container_name /bin/bash 
        else
            # No container with the specified name exists
            echo "Creating and starting a new container $container_name..."
            docker run -it \
                       --gpus all \
                       -e DISPLAY=$DISPLAY \
                       -v $PWD:/workspace/src/$SF_NAME \
                       -v /tmp/.X11-unix/:/tmp/.X11-unix/ \
					   --device=/dev/dri:/dev/dri \
	                   --name $container_name $image_name /bin/bash 
						
        fi
    fi

cd -
