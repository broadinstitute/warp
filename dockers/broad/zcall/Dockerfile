FROM alpine:3.8

ENV TERM=xterm-256color

WORKDIR /usr/gitc

# Install dependencies.
RUN apk --update add bash python2

COPY . .
