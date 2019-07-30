FROM ubuntu:16.04

CMD ["/bin/bash"]

RUN apt-get update

COPY ../refseq_correct.py /usr/refseq_correct.py

COPY ../rules_filter_v3_MuTect.py /usr/rules_filter_v3_MuTect.py

COPY ../rules_filter.py /usr/rules_filter.py

RUN apt-get update && apt-get install -y autoconf automake bison build-essential fail2ban gawk git-core groff groff-base erlang-dev libasound2-dev libavcodec-dev libavutil-dev libavformat-dev libav-tools libavresample-dev libswscale-dev liba52-0.7.4-dev libssl-dev libdb-dev libexpat1-dev libcurl4-openssl-dev libgdbm-dev libgnutls-dev libjpeg-dev libmp3lame-dev libncurses5 libncurses5-dev libperl-dev libogg-dev libsnmp-dev libtiff5-dev libtool libvorbis-dev libx11-dev libzrtpcpp-dev make portaudio19-dev python-dev snmp snmpd subversion unixodbc-dev uuid-dev zlib1g-dev libsqlite3-dev libpcre3-dev libspeex-dev libspeexdsp-dev libldns-dev libedit-dev libladspa-ocaml-dev libmemcached-dev libmp4v2-dev libpq-dev libvlc-dev libv8-dev liblua5.2-dev libyaml-dev libpython-dev odbc-postgresql sendmail unixodbc wget yasm libldap2-dev

ENTRYPOINT ["/bin/bash"]