#C declaration
#       extern vector_t biot_savart(const vector_t v, const coil_t *coil)#
#Arguments
#       *ret    %rdi
#       *coil   %rsi
#       v       8(%rsp)
#Locals
#       B       %ymm0 -> 16(%rsp)
        .section .text
        .p2align 4,,15
        .global biot_savart
        .type   biot_savart, @function
biot_savart:
        subq    $8, %rsp
        vxorpd  %ymm0, %ymm0, %ymm0
        vmovupd %ymm0, 16(%rsp)
        addq    $8, %rsp
        ret
        .size   biot_savart, .-biot_savart
        .section .rodata,"a",@progbits
        .align  8
.pi:
        .double 3.14159265359
