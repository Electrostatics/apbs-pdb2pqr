FROM python:3.7-alpine

LABEL "com.github.actions.name" = "Copyright Action"
LABEL "com.github.actions.description" = "Automatically inserts copyright to the configured directories in a repository"
LABEL "com.github.actions.icon" = "edit"
LABEL "com.github.actions.color" = "gray-dark"

LABEL "repository" = "https://github.com/VinnyBabuManjaly/copyright-action"
LABEL "homepage" = "https://github.com/VinnyBabuManjaly/copyright-action"
LABEL maintainer="Vinny Babu Manjaly <vinnybabumanjaly@gmail.com>"
LABEL version="1.0.0"

# Creates work dir
WORKDIR /app

# Copy python script to the container folder app
COPY copyright_insert.py /app/copyright_insert.py

ENTRYPOINT ["python", "/app/copyright_insert.py"]