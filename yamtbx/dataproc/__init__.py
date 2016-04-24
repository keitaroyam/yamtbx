def software_binning(data, binning, beamxy):
    u = l = 0
    b = r = None
    newshape = [data.shape[0]//binning, data.shape[1]//binning]
    if data.shape[0]%binning > 0:
        res = data.shape[0]%binning
        u = res//2 # upper
        b = -(res - u) # bottom
    if data.shape[1]%binning > 0:
        res = data.shape[1]%binning
        l = res//2
        r = -(res - l)

    newdata = data[u:b, l:r]

    # Reference: http://stackoverflow.com/a/8090605
    newdata = newdata.reshape(newshape[0], binning, newshape[1], -1).sum(3).sum(1).astype(data.dtype)

    # Beam center shift
    newbeamxy = (beamxy[0]-u)/binning, (beamxy[1]-l)/binning

    return newdata, newbeamxy
# software_binning()
