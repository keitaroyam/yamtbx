# LIBTBX_SET_DISPATCHER_NAME kamo.single_images_integration

from yamtbx.dataproc.auto.command_line import single_images_integration

if __name__ == "__main__":
    import sys
    single_images_integration.run_from_args(sys.argv[1:])
