STATE2ABBREV = {
    'Alabama': 'AL',
    'Alaska': 'AK',
    'American Samoa': 'AS',
    'Arizona': 'AZ',
    'Arkansas': 'AR',
    'California': 'CA',
    'Colorado': 'CO',
    'Connecticut': 'CT',
    'Delaware': 'DE',
    'District of Columbia': 'DC',
    'Florida': 'FL',
    'Georgia': 'GA',
    'Guam': 'GU',
    'Hawaii': 'HI',
    'Idaho': 'ID',
    'Illinois': 'IL',
    'Indiana': 'IN',
    'Iowa': 'IA',
    'Kansas': 'KS',
    'Kentucky': 'KY',
    'Louisiana': 'LA',
    'Maine': 'ME',
    'Maryland': 'MD',
    'Massachusetts': 'MA',
    'Michigan': 'MI',
    'Minnesota': 'MN',
    'Mississippi': 'MS',
    'Missouri': 'MO',
    'Montana': 'MT',
    'Nebraska': 'NE',
    'Nevada': 'NV',
    'New Hampshire': 'NH',
    'New Jersey': 'NJ',
    'New Mexico': 'NM',
    'New York': 'NY',
    'North Carolina': 'NC',
    'North Dakota': 'ND',
    'Northern Mariana Islands':'MP',
    'Ohio': 'OH',
    'Oklahoma': 'OK',
    'Oregon': 'OR',
    'Pennsylvania': 'PA',
    'Puerto Rico': 'PR',
    'Rhode Island': 'RI',
    'South Carolina': 'SC',
    'South Dakota': 'SD',
    'Tennessee': 'TN',
    'Texas': 'TX',
    'Utah': 'UT',
    'Vermont': 'VT',
    'Virgin Islands': 'VI',
    'Virginia': 'VA',
    'Washington': 'WA',
    'West Virginia': 'WV',
    'Wisconsin': 'WI',
    'Wyoming': 'WY'
}


COUNTY_CORRECTIONS = {
    'Parish': 'County',
    'Manhattan': 'New York',
    'Brooklyn': 'Kings',
    'Staton Island': 'Richmond',
    'New Orleans': 'Orleans',
    'Pittsburgh': 'Alleghany',
    'Boston': 'Suffolk',
    'East Haven': 'New Haven',
}

GENE2POS = {
            '5UTR': {'start': 0, 'end': 265},
            'ORF1a': {'start': 265, 'end': 13466},
            'ORF1b': {'start': 13467, 'end': 21555},
            'S': {'start': 21562, 'end': 25384},
            'ORF3a': {'start': 25392, 'end': 26220},
            'E': {'start': 26244, 'end': 26472},
            'M': {'start': 26522, 'end': 27191},
            'ORF6': {'start': 27201, 'end': 27387},
            'ORF7a': {'start': 27393, 'end': 27759},
            'ORF7b': {'start': 27755, 'end': 27887},
            'ORF8': {'start': 27893, 'end': 28259},
            'N': {'start': 28273, 'end': 29533},
            'ORF10': {'start': 29557, 'end': 29674},
            '3UTR': {'start': 29674, 'end': 29902}
           }