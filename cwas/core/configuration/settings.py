"""
Default category domain information for CWAS project.
"""

# This domain dictionary is incomplete.
# 'conservation' domains are added by BigWig file keys.
# 'gene_list' domains are added by gene matrix.
# 'region' domains are added by BED file keys.
default_domains = {
    'variant_type': [
        'All',
        'SNV',
        'Indel',
    ],
    'conservation': [
        'All',
    ],
    'gene_list': [
        'Any',
    ],
    'gencode': [  # GENCODE annotation categories
        'Any',
        'CodingRegion',
        'FrameshiftRegion',
        'InFrameRegion',
        'SilentRegion',
        'LoFRegion',
        'MissenseHVARDRegionSimple',
        'MissenseRegion',
        'NoncodingRegion',
        'SpliceSiteNoncanonRegion',
        'IntronRegion',
        'PromoterRegion',
        'IntergenicRegion',
        'UTRsRegion',
        'AntisenseRegion',
        'lincRnaRegion',
        'OtherTranscriptRegion',
    ],
    'region': [  # Custom annotation categories
        'Any',
    ]
}

# A category (domain combination) that includes a redundant domain pair will be
# excluded from CWAS analysis.
redundant_domain_pairs = {
    ('variant_type', 'gencode'): {
        ('All', 'FrameshiftRegion'),
        ('All', 'InFrameRegion'),
        ('All', 'MissenseHVARDRegionSimple'),
        ('All', 'MissenseRegion'),
        ('All', 'SilentRegion'),
        ('Indel', 'MissenseHVARDRegionSimple'),
        ('Indel', 'MissenseRegion'),
        ('Indel', 'SilentRegion'),
        ('SNV', 'FrameshiftRegion'),
        ('SNV', 'InFrameRegion'),
    },
    ('gene_list', 'gencode'): {
        ('Any', 'AntisenseRegion'),
        ('Any', 'CodingRegion'),
        ('Any', 'FrameshiftRegion'),
        ('Any', 'InFrameRegion'),
        ('Any', 'LoFRegion'),
        ('Any', 'MissenseHVARDRegionSimple'),
        ('Any', 'MissenseRegion'),
        ('Any', 'SilentRegion'),
        ('Any', 'lincRnaRegion'),
    }
}
