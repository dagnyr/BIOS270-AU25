FROM dagnyr/bioinfo_example:amd64

RUN pip install parasail

# Make sure HOME/bin exists and is writable
RUN mkdir -p $HOME/bin && \
    curl -L https://github.com/rcedgar/reseek/releases/download/v2.7/reseek-v2.7-linux-x86 \
        -o $HOME/bin/reseek && \
    chmod +x $HOME/bin/reseek

# Add to PATH
ENV PATH="$HOME/bin:${PATH}"

