import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional

import requests
from latch.executions import rename_current_execution, report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

sys.stdout.reconfigure(line_buffering=True)

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize(run_name: str) -> str:
    rename_current_execution(str(run_name))

    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_expiration_hours": 0,
            "version": 2,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


@dataclass
class SampleSheet:
    group: str
    replicate: int
    fastq_1: LatchFile
    fastq_2: LatchFile
    control: str


class Reference(Enum):
    hg38 = "GRCh38 (Homo Sapiens hg38)"
    hg19 = "GRCh37 (Homo Sapiens hg19)"


input_construct_samplesheet = metadata._nextflow_metadata.parameters["input"].samplesheet_constructor


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(
    pvc_name: str,
    input: List[SampleSheet],
    outdir: LatchOutputDir,
    multiqc_title: Optional[str],
    save_reference: bool,
    save_merged_fastq: bool,
    save_trimmed: bool,
    save_spikein_aligned: bool,
    save_unaligned: bool,
    save_align_intermed: bool,
    email: Optional[str],
    dump_scale_factors: bool,
    genome_source: str,
    genome: Optional[str],
    latch_genome: Optional[Reference],
    bowtie2: Optional[LatchDir],
    gtf: Optional[LatchFile],
    gene_bed: Optional[LatchFile],
    blacklist: Optional[LatchFile],
    spikein_bowtie2: Optional[LatchFile],
    spikein_fasta: Optional[LatchFile],
    fasta: Optional[LatchFile],
    only_input: bool,
    only_genome: bool,
    only_preqc: bool,
    only_alignment: bool,
    only_filtering: bool,
    only_peak_calling: bool,
    skip_fastqc: bool,
    skip_trimming: bool,
    skip_removeduplicates: bool,
    skip_reporting: bool,
    skip_preseq: bool,
    skip_igv: bool,
    skip_dt_qc: bool,
    skip_heatmaps: bool,
    skip_peak_qc: bool,
    skip_multiqc: bool,
    remove_mitochondrial_reads: bool,
    run_name: str,
    mito_name: Optional[str],
    dedup_target_reads: bool,
    remove_linear_duplicates: bool,
    macs2_pvalue: Optional[float],
    singularity_pull_docker_container: bool,
    multiqc_methods_description: Optional[str],
    spikein_genome: Optional[str],
    clip_r1: Optional[int],
    clip_r2: Optional[int],
    three_prime_clip_r1: Optional[int],
    three_prime_clip_r2: Optional[int],
    trim_nextseq: Optional[int],
    minimum_alignment_q_score: Optional[int],
    end_to_end: bool,
    normalisation_mode: Optional[str],
    normalisation_binsize: Optional[int],
    peakcaller: Optional[str],
    use_control: bool,
    extend_fragments: bool,
    igg_scale_factor: Optional[float],
    seacr_peak_threshold: Optional[float],
    seacr_norm: Optional[str],
    seacr_stringent: Optional[str],
    macs2_qvalue: Optional[float],
    macs_gsize: Optional[float],
    macs2_narrow_peak: bool,
    macs2_broad_cutoff: Optional[float],
    consensus_peak_mode: Optional[str],
    replicate_threshold: Optional[float],
    igv_show_gene_names: bool,
    dt_qc_bam_binsize: Optional[int],
    dt_qc_corr_method: Optional[str],
    dt_heatmap_gene_beforelen: Optional[int],
    dt_heatmap_gene_bodylen: Optional[int],
    dt_heatmap_gene_afterlen: Optional[int],
    dt_heatmap_peak_beforelen: Optional[int],
    dt_heatmap_peak_afterlen: Optional[int],
    dt_calc_all_matrix: bool,
    min_frip_overlap: Optional[float],
    min_peak_overlap: Optional[float],
    igv_sort_by_groups: bool,
) -> None:
    shared_dir = Path("/nf-workdir")

    input_samplesheet = input_construct_samplesheet(input)

    ignore_list = [
        "latch",
        ".latch",
        ".git",
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

    if genome_source == "latch_genome_source":
        fasta = os.path.join(
            "s3://latch-public/test-data/35729/cutandrun/",
            latch_genome.name,
            latch_genome.name + ".fa",
        )
        gtf = os.path.join(
            "s3://latch-public/test-data/35729/cutandrun/",
            latch_genome.name,
            latch_genome.name + ".refGene.gtf",
        )
        bowtie2 = os.path.join("s3://latch-public/test-data/35729/cutandrun/", latch_genome.name, "bowtie2")

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
        "-resume",
        *get_flag("input", input_samplesheet),
        *get_flag("outdir", LatchOutputDir(outdir.remote_path + "/" + run_name)),
        *get_flag("multiqc_title", multiqc_title),
        *get_flag("save_reference", save_reference),
        *get_flag("save_merged_fastq", save_merged_fastq),
        *get_flag("save_trimmed", save_trimmed),
        *get_flag("save_spikein_aligned", save_spikein_aligned),
        *get_flag("save_unaligned", save_unaligned),
        *get_flag("save_align_intermed", save_align_intermed),
        *get_flag("email", email),
        *get_flag("dump_scale_factors", dump_scale_factors),
        *get_flag("genome", genome),
        *get_flag("bowtie2", bowtie2),
        *get_flag("gtf", gtf),
        *get_flag("gene_bed", gene_bed),
        *get_flag("blacklist", blacklist),
        *get_flag("spikein_genome", spikein_genome),
        *get_flag("spikein_bowtie2", spikein_bowtie2),
        *get_flag("spikein_fasta", spikein_fasta),
        *get_flag("fasta", fasta),
        *get_flag("only_input", only_input),
        *get_flag("only_genome", only_genome),
        *get_flag("only_preqc", only_preqc),
        *get_flag("only_alignment", only_alignment),
        *get_flag("only_filtering", only_filtering),
        *get_flag("only_peak_calling", only_peak_calling),
        *get_flag("skip_fastqc", skip_fastqc),
        *get_flag("skip_trimming", skip_trimming),
        *get_flag("skip_removeduplicates", skip_removeduplicates),
        *get_flag("skip_reporting", skip_reporting),
        *get_flag("skip_preseq", skip_preseq),
        *get_flag("skip_igv", skip_igv),
        *get_flag("skip_dt_qc", skip_dt_qc),
        *get_flag("skip_heatmaps", skip_heatmaps),
        *get_flag("skip_peak_qc", skip_peak_qc),
        *get_flag("skip_multiqc", skip_multiqc),
        *get_flag("clip_r1", clip_r1),
        *get_flag("clip_r2", clip_r2),
        *get_flag("three_prime_clip_r1", three_prime_clip_r1),
        *get_flag("three_prime_clip_r2", three_prime_clip_r2),
        *get_flag("trim_nextseq", trim_nextseq),
        *get_flag("minimum_alignment_q_score", minimum_alignment_q_score),
        *get_flag("remove_mitochondrial_reads", remove_mitochondrial_reads),
        *get_flag("mito_name", mito_name),
        *get_flag("dedup_target_reads", dedup_target_reads),
        *get_flag("remove_linear_duplicates", remove_linear_duplicates),
        *get_flag("end_to_end", end_to_end),
        *get_flag("normalisation_mode", normalisation_mode),
        *get_flag("normalisation_binsize", normalisation_binsize),
        *get_flag("peakcaller", peakcaller),
        *get_flag("use_control", use_control),
        *get_flag("extend_fragments", extend_fragments),
        *get_flag("igg_scale_factor", igg_scale_factor),
        *get_flag("seacr_peak_threshold", seacr_peak_threshold),
        *get_flag("seacr_norm", seacr_norm),
        *get_flag("seacr_stringent", seacr_stringent),
        *get_flag("macs2_qvalue", macs2_qvalue),
        *get_flag("macs2_pvalue", macs2_pvalue),
        *get_flag("macs_gsize", macs_gsize),
        *get_flag("macs2_narrow_peak", macs2_narrow_peak),
        *get_flag("macs2_broad_cutoff", macs2_broad_cutoff),
        *get_flag("consensus_peak_mode", consensus_peak_mode),
        *get_flag("replicate_threshold", replicate_threshold),
        *get_flag("igv_show_gene_names", igv_show_gene_names),
        *get_flag("dt_qc_bam_binsize", dt_qc_bam_binsize),
        *get_flag("dt_qc_corr_method", dt_qc_corr_method),
        *get_flag("dt_heatmap_gene_beforelen", dt_heatmap_gene_beforelen),
        *get_flag("dt_heatmap_gene_bodylen", dt_heatmap_gene_bodylen),
        *get_flag("dt_heatmap_gene_afterlen", dt_heatmap_gene_afterlen),
        *get_flag("dt_heatmap_peak_beforelen", dt_heatmap_peak_beforelen),
        *get_flag("dt_heatmap_peak_afterlen", dt_heatmap_peak_afterlen),
        *get_flag("dt_calc_all_matrix", dt_calc_all_matrix),
        *get_flag("min_frip_overlap", min_frip_overlap),
        *get_flag("min_peak_overlap", min_peak_overlap),
        *get_flag("igv_sort_by_groups", igv_sort_by_groups),
        *get_flag("singularity_pull_docker_container", singularity_pull_docker_container),
        *get_flag("multiqc_methods_description", multiqc_methods_description),
    ]

    print("Launching Nextflow Runtime")
    print(" ".join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
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

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ["du", "-sb", str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60,
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print("Failed to compute storage size: Operation timed out after 5 minutes.")
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)
