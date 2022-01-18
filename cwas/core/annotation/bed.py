from collections import deque

import pysam


# TODO: Make this code much clearer
def annotate(in_vcf_gz_path: str, out_vcf_path: str, annot_bed_path: str):
    chroms = [f"chr{n}" for n in range(1, 23)]

    with pysam.TabixFile(in_vcf_gz_path) as in_vcf_file, pysam.TabixFile(
        annot_bed_path
    ) as annot_bed_file, open(out_vcf_path, "w") as out_vcf_file:
        # Make and write headers
        vcf_headers = in_vcf_file.header
        annot_key_str = annot_bed_file.header[0].split("=")[1]
        annot_info_header = f"##INFO=<ID=ANNOT,Key={annot_key_str}>"
        vcf_headers.append(annot_info_header)
        vcf_headers[-1], vcf_headers[-2] = (
            vcf_headers[-2],
            vcf_headers[-1],
        )  # Swap

        for vcf_header in vcf_headers:
            print(vcf_header, file=out_vcf_file)

        # Annotate by the input BED file
        for chrom in chroms:
            var_iter = in_vcf_file.fetch(chrom, parser=pysam.asTuple())
            bed_iter = annot_bed_file.fetch(chrom, parser=pysam.asTuple())
            bed_memory = deque()
            variant = next(var_iter, None)

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

                # 2. Continuously iterate over the list BED coordinates and check
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

                print(str(variant) + f";ANNOT={annot_int}", file=out_vcf_file)
                variant = next(var_iter, None)
