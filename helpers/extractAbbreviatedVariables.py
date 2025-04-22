import netCDF4 as nc
toexclude = ['XLAT', 'XLONG', 'PH', 'PHB', 'P', 'PB', 'W', 'QRAIN']
toexclude = ['W']

with nc.Dataset("/Users/nbarron/work/AR-20211024/600m/wrfout_d01_002990.nc") as src, nc.Dataset("/Users/nbarron/work/AR-20211024/600m/abbreviated/wrfout_d01_002990.nc", "w") as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if name in toexclude:
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = src[name][:]
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)