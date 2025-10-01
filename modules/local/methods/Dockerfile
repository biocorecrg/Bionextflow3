FROM python:3.12.11

ARG user=node

RUN pip install pyyaml==6.0.2
RUN apt-get update; \
    apt-get full-upgrade -y; \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
    ; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*


# Add non-root user
RUN useradd -ms /bin/bash ${user}
USER ${user}

# we need a .bashrc  for install.sh to work
RUN touch ~/.bashrc && chmod +x ~/.bashrc

ENV NVM_VERSION v0.40.3
ENV NODE_VERSION v22.20.0

RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/${NVM_VERSION}/install.sh | bash ;

# nvm
RUN echo 'export NVM_DIR="$HOME/.nvm"'                                       >> "$HOME/.bashrc"
RUN echo '[ -s "$NVM_DIR/nvm.sh" ] && . "$NVM_DIR/nvm.sh"  # This loads nvm' >> "$HOME/.bashrc"
RUN echo '[ -s "$NVM_DIR/bash_completion" ] && . "$NVM_DIR/bash_completion" # This loads nvm bash_completion' >> "$HOME/.bashrc"


# nodejs and tools
RUN bash -c 'source $HOME/.nvm/nvm.sh   && \
    nvm install $NODE_VERSION && nvm use $NODE_VERSION && nvm alias default $NODE_VERSION && \
    npm install --global @citation-js/core@^0.7.18 \
     @citation-js/plugin-csl@^0.7.18 \
     @citation-js/plugin-doi@^0.7.18 \
     turndown@^7.2.1 \
     js-yaml && \
    NPM_GLOBAL_BIN="$HOME/.nvm/versions/node/$NODE_VERSION/bin" && \
    echo "export PATH=$NPM_GLOBAL_BIN:\$PATH" >> $HOME/.bashrc && \
    echo "export PATH=$NPM_GLOBAL_BIN:\$PATH" >> $HOME/.profile';

ENV PATH=/home/node/.nvm/versions/node/$NODE_VERSION/bin:$PATH
ENV NODE_PATH=/home/node/.nvm/versions/node/$NODE_VERSION/lib/node_modules

# then  the rest of the dockerfile
WORKDIR /home/${user}
