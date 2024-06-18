from dataclasses import dataclass
from enum import Enum
import os
import subprocess
import requests
import shutil
from pathlib import Path
import typing
import typing_extensions

from latch.resources.workflow import workflow
from latch.resources.tasks import nextflow_runtime_task, custom_task
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.ldata.path import LPath
from latch_cli.nextflow.workflow import get_flag
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.utils import urljoins
from latch.types import metadata
from flytekit.core.annotation import FlyteAnnotation

from latch_cli.services.register.utils import import_module_by_path

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata

@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_gib": 100,
        }
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]






@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(pvc_name: str, input: LatchFile, outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], multiqc_title: typing.Optional[str], save_reference: typing.Optional[bool], save_merged_fastq: typing.Optional[bool], save_trimmed: typing.Optional[bool], save_spikein_aligned: typing.Optional[bool], save_unaligned: typing.Optional[bool], save_align_intermed: typing.Optional[bool], email: typing.Optional[str], dump_scale_factors: typing.Optional[bool], genome: typing.Optional[str], bowtie2: typing.Optional[str], gtf: typing.Optional[str], gene_bed: typing.Optional[str], blacklist: typing.Optional[str], spikein_bowtie2: typing.Optional[str], spikein_fasta: typing.Optional[str], fasta: typing.Optional[LatchFile], only_input: typing.Optional[bool], only_genome: typing.Optional[bool], only_preqc: typing.Optional[bool], only_alignment: typing.Optional[bool], only_filtering: typing.Optional[bool], only_peak_calling: typing.Optional[bool], skip_fastqc: typing.Optional[bool], skip_trimming: typing.Optional[bool], skip_removeduplicates: typing.Optional[bool], skip_reporting: typing.Optional[bool], skip_preseq: typing.Optional[bool], skip_igv: typing.Optional[bool], skip_dt_qc: typing.Optional[bool], skip_heatmaps: typing.Optional[bool], skip_peak_qc: typing.Optional[bool], skip_multiqc: typing.Optional[bool], remove_mitochondrial_reads: typing.Optional[bool], mito_name: typing.Optional[str], dedup_target_reads: typing.Optional[bool], remove_linear_duplicates: typing.Optional[bool], macs2_pvalue: typing.Optional[float], singularity_pull_docker_container: typing.Optional[bool], multiqc_methods_description: typing.Optional[str], spikein_genome: typing.Optional[str], clip_r1: typing.Optional[int], clip_r2: typing.Optional[int], three_prime_clip_r1: typing.Optional[int], three_prime_clip_r2: typing.Optional[int], trim_nextseq: typing.Optional[int], minimum_alignment_q_score: typing.Optional[int], end_to_end: typing.Optional[bool], normalisation_mode: typing.Optional[str], normalisation_binsize: typing.Optional[int], peakcaller: typing.Optional[str], use_control: typing.Optional[bool], extend_fragments: typing.Optional[bool], igg_scale_factor: typing.Optional[float], seacr_peak_threshold: typing.Optional[float], seacr_norm: typing.Optional[str], seacr_stringent: typing.Optional[str], macs2_qvalue: typing.Optional[float], macs_gsize: typing.Optional[float], macs2_narrow_peak: typing.Optional[bool], macs2_broad_cutoff: typing.Optional[float], consensus_peak_mode: typing.Optional[str], replicate_threshold: typing.Optional[float], igv_show_gene_names: typing.Optional[bool], dt_qc_bam_binsize: typing.Optional[int], dt_qc_corr_method: typing.Optional[str], dt_heatmap_gene_beforelen: typing.Optional[int], dt_heatmap_gene_bodylen: typing.Optional[int], dt_heatmap_gene_afterlen: typing.Optional[int], dt_heatmap_peak_beforelen: typing.Optional[int], dt_heatmap_peak_afterlen: typing.Optional[int], dt_calc_all_matrix: typing.Optional[bool], min_frip_overlap: typing.Optional[float], min_peak_overlap: typing.Optional[float], igv_sort_by_groups: typing.Optional[bool]) -> None:
    try:
        shared_dir = Path("/nf-workdir")



        ignore_list = [
            "latch",
            ".latch",
            "nextflow",
            ".nextflow",
            "work",
            "results",
            "miniconda",
            "anaconda3",
            "mambaforge",
        ]

        shutil.copytree(
            Path("/root"),
            shared_dir,
            ignore=lambda src, names: ignore_list,
            ignore_dangling_symlinks=True,
            dirs_exist_ok=True,
        )

        cmd = [
            "/root/nextflow",
            "run",
            str(shared_dir / "main.nf"),
            "-work-dir",
            str(shared_dir),
            "-profile",
            "docker",
            "-c",
            "latch.config",
                *get_flag('input', input),
                *get_flag('outdir', outdir),
                *get_flag('multiqc_title', multiqc_title),
                *get_flag('save_reference', save_reference),
                *get_flag('save_merged_fastq', save_merged_fastq),
                *get_flag('save_trimmed', save_trimmed),
                *get_flag('save_spikein_aligned', save_spikein_aligned),
                *get_flag('save_unaligned', save_unaligned),
                *get_flag('save_align_intermed', save_align_intermed),
                *get_flag('email', email),
                *get_flag('dump_scale_factors', dump_scale_factors),
                *get_flag('genome', genome),
                *get_flag('bowtie2', bowtie2),
                *get_flag('gtf', gtf),
                *get_flag('gene_bed', gene_bed),
                *get_flag('blacklist', blacklist),
                *get_flag('spikein_genome', spikein_genome),
                *get_flag('spikein_bowtie2', spikein_bowtie2),
                *get_flag('spikein_fasta', spikein_fasta),
                *get_flag('fasta', fasta),
                *get_flag('only_input', only_input),
                *get_flag('only_genome', only_genome),
                *get_flag('only_preqc', only_preqc),
                *get_flag('only_alignment', only_alignment),
                *get_flag('only_filtering', only_filtering),
                *get_flag('only_peak_calling', only_peak_calling),
                *get_flag('skip_fastqc', skip_fastqc),
                *get_flag('skip_trimming', skip_trimming),
                *get_flag('skip_removeduplicates', skip_removeduplicates),
                *get_flag('skip_reporting', skip_reporting),
                *get_flag('skip_preseq', skip_preseq),
                *get_flag('skip_igv', skip_igv),
                *get_flag('skip_dt_qc', skip_dt_qc),
                *get_flag('skip_heatmaps', skip_heatmaps),
                *get_flag('skip_peak_qc', skip_peak_qc),
                *get_flag('skip_multiqc', skip_multiqc),
                *get_flag('clip_r1', clip_r1),
                *get_flag('clip_r2', clip_r2),
                *get_flag('three_prime_clip_r1', three_prime_clip_r1),
                *get_flag('three_prime_clip_r2', three_prime_clip_r2),
                *get_flag('trim_nextseq', trim_nextseq),
                *get_flag('minimum_alignment_q_score', minimum_alignment_q_score),
                *get_flag('remove_mitochondrial_reads', remove_mitochondrial_reads),
                *get_flag('mito_name', mito_name),
                *get_flag('dedup_target_reads', dedup_target_reads),
                *get_flag('remove_linear_duplicates', remove_linear_duplicates),
                *get_flag('end_to_end', end_to_end),
                *get_flag('normalisation_mode', normalisation_mode),
                *get_flag('normalisation_binsize', normalisation_binsize),
                *get_flag('peakcaller', peakcaller),
                *get_flag('use_control', use_control),
                *get_flag('extend_fragments', extend_fragments),
                *get_flag('igg_scale_factor', igg_scale_factor),
                *get_flag('seacr_peak_threshold', seacr_peak_threshold),
                *get_flag('seacr_norm', seacr_norm),
                *get_flag('seacr_stringent', seacr_stringent),
                *get_flag('macs2_qvalue', macs2_qvalue),
                *get_flag('macs2_pvalue', macs2_pvalue),
                *get_flag('macs_gsize', macs_gsize),
                *get_flag('macs2_narrow_peak', macs2_narrow_peak),
                *get_flag('macs2_broad_cutoff', macs2_broad_cutoff),
                *get_flag('consensus_peak_mode', consensus_peak_mode),
                *get_flag('replicate_threshold', replicate_threshold),
                *get_flag('igv_show_gene_names', igv_show_gene_names),
                *get_flag('dt_qc_bam_binsize', dt_qc_bam_binsize),
                *get_flag('dt_qc_corr_method', dt_qc_corr_method),
                *get_flag('dt_heatmap_gene_beforelen', dt_heatmap_gene_beforelen),
                *get_flag('dt_heatmap_gene_bodylen', dt_heatmap_gene_bodylen),
                *get_flag('dt_heatmap_gene_afterlen', dt_heatmap_gene_afterlen),
                *get_flag('dt_heatmap_peak_beforelen', dt_heatmap_peak_beforelen),
                *get_flag('dt_heatmap_peak_afterlen', dt_heatmap_peak_afterlen),
                *get_flag('dt_calc_all_matrix', dt_calc_all_matrix),
                *get_flag('min_frip_overlap', min_frip_overlap),
                *get_flag('min_peak_overlap', min_peak_overlap),
                *get_flag('igv_sort_by_groups', igv_sort_by_groups),
                *get_flag('singularity_pull_docker_container', singularity_pull_docker_container),
                *get_flag('multiqc_methods_description', multiqc_methods_description)
        ]

        print("Launching Nextflow Runtime")
        print(' '.join(cmd))
        print(flush=True)

        env = {
            **os.environ,
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms2048M -Xmx8G -XX:ActiveProcessorCount=4",
            "K8S_STORAGE_CLAIM_NAME": pvc_name,
            "NXF_DISABLE_CHECK_LATEST": "true",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(urljoins("latch:///your_log_dir/nf_nf_core_cutandrun", name, "nextflow.log"))
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)



@workflow(metadata._nextflow_metadata)
def nf_nf_core_cutandrun(input: LatchFile, outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], multiqc_title: typing.Optional[str], save_reference: typing.Optional[bool], save_merged_fastq: typing.Optional[bool], save_trimmed: typing.Optional[bool], save_spikein_aligned: typing.Optional[bool], save_unaligned: typing.Optional[bool], save_align_intermed: typing.Optional[bool], email: typing.Optional[str], dump_scale_factors: typing.Optional[bool], genome: typing.Optional[str], bowtie2: typing.Optional[str], gtf: typing.Optional[str], gene_bed: typing.Optional[str], blacklist: typing.Optional[str], spikein_bowtie2: typing.Optional[str], spikein_fasta: typing.Optional[str], fasta: typing.Optional[LatchFile], only_input: typing.Optional[bool], only_genome: typing.Optional[bool], only_preqc: typing.Optional[bool], only_alignment: typing.Optional[bool], only_filtering: typing.Optional[bool], only_peak_calling: typing.Optional[bool], skip_fastqc: typing.Optional[bool], skip_trimming: typing.Optional[bool], skip_removeduplicates: typing.Optional[bool], skip_reporting: typing.Optional[bool], skip_preseq: typing.Optional[bool], skip_igv: typing.Optional[bool], skip_dt_qc: typing.Optional[bool], skip_heatmaps: typing.Optional[bool], skip_peak_qc: typing.Optional[bool], skip_multiqc: typing.Optional[bool], remove_mitochondrial_reads: typing.Optional[bool], mito_name: typing.Optional[str], dedup_target_reads: typing.Optional[bool], remove_linear_duplicates: typing.Optional[bool], macs2_pvalue: typing.Optional[float], singularity_pull_docker_container: typing.Optional[bool], multiqc_methods_description: typing.Optional[str], spikein_genome: typing.Optional[str] = 'K12-MG1655', clip_r1: typing.Optional[int] = 0, clip_r2: typing.Optional[int] = 0, three_prime_clip_r1: typing.Optional[int] = 0, three_prime_clip_r2: typing.Optional[int] = 0, trim_nextseq: typing.Optional[int] = 0, minimum_alignment_q_score: typing.Optional[int] = 20, end_to_end: typing.Optional[bool] = True, normalisation_mode: typing.Optional[str] = 'Spikein', normalisation_binsize: typing.Optional[int] = 50, peakcaller: typing.Optional[str] = 'seacr', use_control: typing.Optional[bool] = True, extend_fragments: typing.Optional[bool] = True, igg_scale_factor: typing.Optional[float] = 0.5, seacr_peak_threshold: typing.Optional[float] = 0.05, seacr_norm: typing.Optional[str] = 'non', seacr_stringent: typing.Optional[str] = 'stringent', macs2_qvalue: typing.Optional[float] = 0.01, macs_gsize: typing.Optional[float] = 2700000000, macs2_narrow_peak: typing.Optional[bool] = True, macs2_broad_cutoff: typing.Optional[float] = 0.1, consensus_peak_mode: typing.Optional[str] = 'group', replicate_threshold: typing.Optional[float] = 1, igv_show_gene_names: typing.Optional[bool] = True, dt_qc_bam_binsize: typing.Optional[int] = 500, dt_qc_corr_method: typing.Optional[str] = 'pearson', dt_heatmap_gene_beforelen: typing.Optional[int] = 3000, dt_heatmap_gene_bodylen: typing.Optional[int] = 5000, dt_heatmap_gene_afterlen: typing.Optional[int] = 3000, dt_heatmap_peak_beforelen: typing.Optional[int] = 3000, dt_heatmap_peak_afterlen: typing.Optional[int] = 3000, dt_calc_all_matrix: typing.Optional[bool] = True, min_frip_overlap: typing.Optional[float] = 0.2, min_peak_overlap: typing.Optional[float] = 0.2, igv_sort_by_groups: typing.Optional[bool] = True) -> None:
    """
    nf-core/cutandrun

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(pvc_name=pvc_name, input=input, outdir=outdir, multiqc_title=multiqc_title, save_reference=save_reference, save_merged_fastq=save_merged_fastq, save_trimmed=save_trimmed, save_spikein_aligned=save_spikein_aligned, save_unaligned=save_unaligned, save_align_intermed=save_align_intermed, email=email, dump_scale_factors=dump_scale_factors, genome=genome, bowtie2=bowtie2, gtf=gtf, gene_bed=gene_bed, blacklist=blacklist, spikein_genome=spikein_genome, spikein_bowtie2=spikein_bowtie2, spikein_fasta=spikein_fasta, fasta=fasta, only_input=only_input, only_genome=only_genome, only_preqc=only_preqc, only_alignment=only_alignment, only_filtering=only_filtering, only_peak_calling=only_peak_calling, skip_fastqc=skip_fastqc, skip_trimming=skip_trimming, skip_removeduplicates=skip_removeduplicates, skip_reporting=skip_reporting, skip_preseq=skip_preseq, skip_igv=skip_igv, skip_dt_qc=skip_dt_qc, skip_heatmaps=skip_heatmaps, skip_peak_qc=skip_peak_qc, skip_multiqc=skip_multiqc, clip_r1=clip_r1, clip_r2=clip_r2, three_prime_clip_r1=three_prime_clip_r1, three_prime_clip_r2=three_prime_clip_r2, trim_nextseq=trim_nextseq, minimum_alignment_q_score=minimum_alignment_q_score, remove_mitochondrial_reads=remove_mitochondrial_reads, mito_name=mito_name, dedup_target_reads=dedup_target_reads, remove_linear_duplicates=remove_linear_duplicates, end_to_end=end_to_end, normalisation_mode=normalisation_mode, normalisation_binsize=normalisation_binsize, peakcaller=peakcaller, use_control=use_control, extend_fragments=extend_fragments, igg_scale_factor=igg_scale_factor, seacr_peak_threshold=seacr_peak_threshold, seacr_norm=seacr_norm, seacr_stringent=seacr_stringent, macs2_qvalue=macs2_qvalue, macs2_pvalue=macs2_pvalue, macs_gsize=macs_gsize, macs2_narrow_peak=macs2_narrow_peak, macs2_broad_cutoff=macs2_broad_cutoff, consensus_peak_mode=consensus_peak_mode, replicate_threshold=replicate_threshold, igv_show_gene_names=igv_show_gene_names, dt_qc_bam_binsize=dt_qc_bam_binsize, dt_qc_corr_method=dt_qc_corr_method, dt_heatmap_gene_beforelen=dt_heatmap_gene_beforelen, dt_heatmap_gene_bodylen=dt_heatmap_gene_bodylen, dt_heatmap_gene_afterlen=dt_heatmap_gene_afterlen, dt_heatmap_peak_beforelen=dt_heatmap_peak_beforelen, dt_heatmap_peak_afterlen=dt_heatmap_peak_afterlen, dt_calc_all_matrix=dt_calc_all_matrix, min_frip_overlap=min_frip_overlap, min_peak_overlap=min_peak_overlap, igv_sort_by_groups=igv_sort_by_groups, singularity_pull_docker_container=singularity_pull_docker_container, multiqc_methods_description=multiqc_methods_description)

