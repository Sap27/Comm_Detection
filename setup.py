from setuptools import setup, find_packages

setup(
    name='community_detection_1',  # This should be the name of your package
    version='0.1',
    packages=find_packages(include=['community_detection_1', 'community_detection_1.*']),
    install_requires=[
        'numpy',
        'networkx',
        'pandas',
        'matplotlib',
        'infomap',
        'scikit-learn',
        'argparse'
    ],
    
)
entry_points={
    'console_scripts': [
        'community-detect = community_detection_1.main:main',
    ],
}

"""entry_points={
        'console_scripts': [
            'community_detection_1=community_detection_1.community_detection_main:main',  # This maps to your main function
        ],
    },"""