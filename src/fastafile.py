"""
Class that wraps a file descriptor for a fasta file to read genome sequences easily.
"""
import os
import gzip


class FastaFile:
    def __init__(self, fasta_file_path: str):
        """
        :param fasta_file_path: This file path should end with '.fa' or '.fa.gz'.
        """
        if not fasta_file_path.endswith('.fa.gz') and not fasta_file_path.endswith('fa'):
            raise ValueError(f'Invalid file format. The file extension should be .fa or .fa.gz.')

        if fasta_file_path.endswith('.gz'):
            is_gzip = True
        else:
            is_gzip = False

        self._fasta_file = gzip.open(fasta_file_path, 'rt') if is_gzip else open(fasta_file_path, 'r')
        self._idx_info_dict = {}  # Key: Chrom ID, Value: Dictionary contains fields of faidx

        # Parse a fasta index if available
        fasta_idx_path = fasta_file_path.replace('.gz', '.fai') if is_gzip else fasta_file_path + '.fai'

        if not os.path.isfile(fasta_idx_path):
            raise FileNotFoundError('A fai index file of the input fasta file cannot be found. '
                                    'Make the index file.')

        with open(fasta_idx_path) as fasta_idx_file:
            for line in fasta_idx_file:
                fields = line.rstrip('\n').split('\t')
                chrom = fields[0]
                size = int(fields[1])
                start_idx = int(fields[2])
                seq_len = int(fields[3])  # Length of a sequence in one line
                line_len = int(fields[4])  # Length of one line with a newline character
                self._idx_info_dict[chrom] = {
                    'size': size,
                    'start_idx': start_idx,
                    'seq_len': seq_len,
                    'line_len': line_len,
                }

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self._fasta_file.close()

    def get_base(self, chrom: str, pos: int) -> str:
        """ Return a base at the input position (0-based) in the input chromosome
        """
        if self._idx_info_dict.get(chrom) is None:
            raise ValueError(f'Invalid chromosome ID "{chrom}"')

        if pos < 0 or pos >= self._idx_info_dict[chrom]['size']:
            raise ValueError('The input position should be in the range [0, chromosome size).')

        base_idx = self._conv_pos_to_idx(chrom, pos)
        self._fasta_file.seek(base_idx)
        base = self._fasta_file.read(1)

        return base

    def get_seq(self, chrom: str, start: int, end: int) -> str:
        """ Return a sequence at the input position (start, end; 0-based) in the input chromosome
        """
        if self._idx_info_dict.get(chrom) is None:
            raise ValueError(f'Invalid chromosome ID "{chrom}"')

        if start >= end:
            raise ValueError('The end position must be larger than the start position.')

        if start < 0 or end > self._idx_info_dict[chrom]['size']:
            raise ValueError('The sequence should be in the range [0, chromosome size].')

        start_idx = self._conv_pos_to_idx(chrom, start)
        end_idx = self._conv_pos_to_idx(chrom, end)
        self._fasta_file.seek(start_idx)
        seq = self._fasta_file.read(end_idx - start_idx)
        seq = seq.replace('\n', '')  # Remove new lines

        return seq

    def _conv_pos_to_idx(self, chrom: str, pos: int) -> int:
        """ Convert the input position on the chromosome into the index on the fasta file """
        return self._idx_info_dict[chrom]['start_idx'] + \
               pos // self._idx_info_dict[chrom]['seq_len'] * self._idx_info_dict[chrom]['line_len'] + \
               pos % self._idx_info_dict[chrom]['seq_len']
