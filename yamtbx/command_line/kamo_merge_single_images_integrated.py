# LIBTBX_SET_DISPATCHER_NAME kamo.merge_single_images_integrated

from yamtbx.dataproc.auto.command_line import merge_single_images_integrated

if __name__ == "__main__":
    import sys
    merge_single_images_integrated.run_from_args(sys.argv[1:])
