FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install coverage && \
    pip install pathos

# -----------------------------------------

# make TopHat2 working dir
RUN mkdir /kb/deployment/bin/TopHat2
# -----------------------------------------

WORKDIR /kb/module

# download Bowtie2 (ver 2.3.2) and copy to Tophat2 dir
RUN VERSION='2.3.2' \
    && wget --no-verbose "http://sourceforge.net/projects/bowtie-bio/files/bowtie2/${VERSION}/bowtie2-${VERSION}-source.zip" \
    && unzip -q bowtie2-${VERSION}-source.zip \
    && cd bowtie2-${VERSION} \
    && make NO_TBB=1 \
    && cp bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s \
          bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s /kb/deployment/bin/TopHat2 \
    && cd .. \
    && rm -rf bowtie2-${VERSION}*
# -----------------------------------------

# download TopHat2 (ver 2.1.1) and copy to Tophat2 dir
RUN VERSION='2.1.1' \
    && wget --no-verbose "https://ccb.jhu.edu/software/tophat/downloads/tophat-${VERSION}.Linux_x86_64.tar.gz" \
    && tar -xzvf tophat-${VERSION}.Linux_x86_64.tar.gz \
    && rm tophat-${VERSION}.Linux_x86_64.tar.gz \
    && cd tophat-${VERSION}.Linux_x86_64 \
    && cp bam2fastx bam_merge bed_to_juncs contig_to_chr_coords fix_map_ordering gtf_juncs gtf_to_fasta \
          juncs_db long_spanning_reads map2gtf prep_reads sam_juncs samtools_0.1.18 segment_juncs sra_to_solid \
          tophat tophat2 tophat-fusion-post tophat_reports /kb/deployment/bin/TopHat2 \
    && cd .. \
    && rm -rf tophat-${VERSION}*
# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
