import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requirements = [
"atomicwrites==1.4.0",
"attrs==21.4.0",
"colorama==0.4.4",
"iniconfig==1.1.1",
"nmrglue==0.8",
"numpy==1.21.4",
"packaging==21.3",
"pandas==1.4.0",
"pluggy==1.0.0",
"py==1.11.0",
"pyparsing==3.0.7",
"pytest==6.2.5",
"python-dateutil==2.8.2",
"pytz==2021.3",
"scipy==1.7.3",
"six==1.12.0",
"toml==0.10.2",
"rdkit-pypi>=2021.9.4",
"requests==2.27.1",
"patool==1.12"
]

setuptools.setup(
    name='apvalidation',
    version='0.5.27',
    author='liningtonlabs',
    author_email='liningtonlabstest@gmail.com',
    description='Testing installation of Package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/liningtonlab/ap_validation',
    license='MIT',
    packages=['apvalidation'],
    install_requires=requirements,
)
