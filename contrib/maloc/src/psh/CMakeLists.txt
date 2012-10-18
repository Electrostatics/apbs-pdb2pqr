add_items(
    SOURCES
    psh.c
    vcom.c
    vmp.c
    vmpi.c
)

add_items(
    EXTERNAL_HEADERS
    maloc/psh.h
    maloc/vcom.h
    maloc/vmp.h
    maloc/vmpi.h
)

add_items(
    INTERNAL_HEADERS
    psh_p.h
    vcom_p.h
    vmpi_p.h
    vmp_p.h
)
