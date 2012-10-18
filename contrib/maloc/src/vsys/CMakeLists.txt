add_items(
    SOURCES
    vio.c
    viofb.c
    vmem.c
    vnm.c
    vnmfb.c
    vset.c
    ziofb.c
)

add_items(
    EXTERNAL_HEADERS
    maloc/vio.h
    maloc/vmem.h
    maloc/vnm.h
    maloc/vset.h
    maloc/vsys.h
)

add_items(
    INTERNAL_HEADERS
    vio_p.h
    vmem_p.h
    vnm_p.h
    vset_p.h
)

