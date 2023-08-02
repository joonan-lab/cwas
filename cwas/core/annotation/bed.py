from collections import deque
import multiprocessing as mp ## test multiprocessing
import pysam
import itertools
import gzip

class annotate:
    def __init__(self, in_vcf_gz_path: str,
                 out_vcf_path: str, annot_bed_path: str, num_proc: str) -> None:
        self._in_vcf_gz_path = in_vcf_gz_path
        self._out_vcf_path = out_vcf_path
        self._annot_bed_path = annot_bed_path
        self._num_proc = num_proc if num_proc < 23 else 22

    @property
    def in_vcf_gz_path(self) -> str:
        return self._in_vcf_gz_path
    @property
    def out_vcf_path(self) -> str:
        return self._out_vcf_path
    @property
    def annot_bed_path(self) -> str:
        return self._annot_bed_path
    @property
    def num_proc(self) -> int:
        return self._num_proc

    def bed_custom_annotate(self):
        with pysam.TabixFile(self.in_vcf_gz_path) as in_vcf_file, pysam.TabixFile(
            self.annot_bed_path
        ) as annot_bed_file, gzip.open(self.out_vcf_path, "wt") as out_vcf_file:
            # Make and write headers
            vcf_headers = in_vcf_file.header
            annot_key_str = annot_bed_file.header[0].split("=")[1]
            annot_info_header = f"##INFO=<ID=ANNOT,Key={annot_key_str}>"
            vcf_headers.append(annot_info_header)
            vcf_headers[-1], vcf_headers[-2] = (
                vcf_headers[-2],
                vcf_headers[-1],
            )  # Swap
            
            # Write VCF header
            out_vcf_file.write("\n".join(vcf_headers))

            # Annotate by the input BED file
            vcf_chroms = in_vcf_file.contigs
            self._num_proc = self._num_proc if len(vcf_chroms) > self._num_proc else len(vcf_chroms)

            p = mp.Pool(self._num_proc)

            annot_vcf = p.map(self.chr_annotate, vcf_chroms)
            annot_vcf = itertools.chain.from_iterable(annot_vcf)
            
            p.close()
            p.join()
            
            # Write VCF annotated
            out_vcf_file.write("\n")
            out_vcf_file.write("\n".join(annot_vcf))
            out_vcf_file.close()

    def chr_annotate(self, chrom):
        input_vcf = pysam.TabixFile(self.in_vcf_gz_path)
        annot_bed = pysam.TabixFile(self.annot_bed_path)

        var_iter = input_vcf.fetch(chrom, parser=pysam.asTuple())
        bed_iter = annot_bed.fetch(chrom, parser=pysam.asTuple())
        bed_memory = deque()
        variant = next(var_iter, None) # var_iter의 마지막 -> None
        variant_annot = []

        while variant is not None:
            var_pos = int(variant[1]) - 1  # 1-based -> 0-based
            var_ref = variant[3]
            var_alt = variant[4]

            # Determine a search region for BED coordinates
            if len(var_ref) == 1:  # Insertion or substitution
                region_start = var_pos
                region_end = (
                var_pos + 2 if len(var_alt) > 1 else var_pos + 1
                )
            else:  # Deletion
                region_start = var_pos + 1
                region_end = region_start + len(var_ref) - 1

            # Get an annotation integer
            annot_int = 0
            stop_bed_iter = False

            # 1. Check the memory of previously checked BED coordinates
            while (
                len(bed_memory) > 0
                and int(bed_memory[0][2]) <= region_start
            ):
                # Remove non-overlapped coordinates
                bed_memory.popleft()

            for bed in bed_memory:
                if int(bed[1]) < region_end:  # Overlap
                    annot_int |= int(bed[3])
                else:
                    stop_bed_iter = True
                    break

            # 2. Continuously iterate over the BED coordinates and check
            if not stop_bed_iter:
                bed = next(bed_iter, None)

                while bed is not None:
                    bed_start = int(bed[1])
                    bed_end = int(bed[2])

                    if region_start < bed_end:
                        bed_memory.append(bed)
                        if bed_start < region_end:  # Overlap
                            annot_int |= int(bed[3])
                        else:
                            break

                    bed = next(bed_iter, None)

            variant_annot.append(str(variant) + f";ANNOT={annot_int}")

            #print(str(variant) + f";ANNOT={annot_int}", file=out_vcf)
            variant = next(var_iter, None)

        return variant_annot