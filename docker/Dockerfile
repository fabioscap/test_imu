from giorgiogrisetti/ubuntu20.04_cuda_ros_noetic_srrg2_solver


WORKDIR /workspace/src/

RUN (cd srrg2_solver && git pull && git checkout devel)
    # freeze the repository to a specific commit
    # || git checkout -b freeze130923 219cf67eb323717e85c27117374d6cf39f8054a0
RUN (cd srrg2_core && git pull && git checkout devel)  
RUN (cd srrg2_qgl_viewport && git pull && git checkout devel)
RUN (cd srrg_cmake_modules && git pull)
RUN (cd srrg2_config_visualizer && git pull && git checkout devel)
RUN (cd srrg2_executor && git pull && git checkout devel)

RUN cd /workspace/ && catkin build

# color output
RUN echo 'force_color_prompt=yes' >> /root/.bashrc
RUN sed -i 's/^#force_color_prompt=yes/force_color_prompt=yes/' /root/.bashrc
RUN echo "PS1='${debian_chroot:+($debian_chroot)}\[\033[01;36m\]srrg\[\033[00m\]:\[\033[31m\]\w\[\033[00m\]\$ '" >> /root/.bashrc
RUN echo "export LS_COLORS='di=31:ln=36'" >> /root/.bashrc
RUN echo 'source /workspace/devel/setup.bash' >> /root/.bashrc

# shared folder
ARG SF_NAME="test_imu"
RUN cd /workspace/src/ && mkdir $SF_NAME

# gtsam
RUN cd / && git clone https://github.com/borglab/gtsam && cd gtsam && mkdir build && cd build && cmake ../  -DGTSAM_USE_SYSTEM_EIGEN=ON -DGTSAM_TANGENT_PREINTEGRATION=OFF && make install 

# upgrade matplotlib
RUN pip install pandas
RUN pip install --upgrade matplotlib

WORKDIR /
