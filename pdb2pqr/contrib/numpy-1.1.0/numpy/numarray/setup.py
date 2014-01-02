from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('numarray',parent_package,top_path)

    config.add_data_files('numpy/')

    config.add_extension('_capi',
                         sources=['_capi.c'],
                         )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
