NXF_HOME=/storage1/fs1/dspencer/Active/spencerlab/dnidhi/nf-chromoseq-dragen/launch_dragen/.nextflow \
LSF_DOCKER_VOLUMES="/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active /scratch1/fs1/dspencer:/scratch1/fs1/dspencer" \
bsub -g /dspencer/dnidhi -G compute-dspencer -q general \
-e launch_dragen.err -o launch_dragen.log -n 2 -M 12GB \
-R "select[mem>=12000] span[hosts=1] rusage[mem=12000]" \
-a "docker(mdivr/nextflow)" \
-u dnidhi@wustl.edu \
nextflow run main.nf -profile test,dragen4 \
--test_sheet /storage1/fs1/dspencer/Active/spencerlab/dnidhi/nf-chromoseq-dragen/assets/demux_samplesheet.csv \
--rundir /storage1/fs1/dspencer/Active/spencerlab/dnidhi/chromoseq_batch/121123_LH00195_ChromoSeq/rundir/20231211_LH00195_0065_B22GLYYLT3/ \
-entry DEMUX \
-w /scratch1/fs1/dspencer/dnidhi/work \
-dump-channels \
-resume 4ca39c08-a308-4822-ba81-95ab906603e7