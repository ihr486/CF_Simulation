        .intel_syntax noprefix
        .section .text
        .p2align 4,,15
        .global biot_savart
        .type   biot_savart, @function
biot_savart:
        subq    rsp, 40                 #Align stack to 16bytes
        vzeroall
        vmovapd [rdi], xmm0
        movsd   16[rdi], xmm0           #Destination is zero cleared
        movsd   xmm0, [rsi]
        mulsd   xmm0, .b_scale[rip]
        movsd   32[rsp], xmm0           #32[rsp] = 1E-7*2*pi*coil->R/THETA_DIV
        mov     ecx, dword ptr 16[rsi]  #ecx = coil->N
        mov     eax, ecx
        sub     eax, 1
        cvtsi2sd xmm0, eax
        mulsd   xmm0, .neg_half[rip]
        mulsd   xmm0, 8[rsi]
        movsd   24[rsp], xmm0           #24[rsp] = -0.5 * coil->P * (coil->N - 1)
biot_savart_loop2:
        mov     edx, 8
        lea     r8, .cos_table[rip]
        lea     r9, .sin_table[rip]
biot_savart_loop1:
        vmovapd ymm0, [r8]
        vmovapd ymm1, [r9]
        vmovapd ymm2, 32[r8]
        vmovapd ymm3, 32[r9]
        add     r8, 64
        add     r9, 64

        vbroadcastsd ymm4, 48[rsp]      #ymm4[3:0] = v.x
        movsd   xmm5, 56[rsp]
        vbroadcastsd ymm6, 64[rsp]      #ymm6[3:0] = v.z
        vbroadcastsd ymm7, [rsi]        #ymm7[3:0] = coil->R
        vmovapd ymm10, ymm4
        vmovapd ymm11, ymm6
        subsd   xmm5, 24[rsp]

        vfnmadd231pd ymm4, ymm7, ymm0   #ymm4 = v.x - coil->R * cos_table[3:0]
        vbroadcastsd ymm5, xmm5         #ymm5 = v.y - 24[rsp]
        vfnmadd231pd ymm6, ymm7, ymm1   #ymm6 = v.z - coil->R * sin_table[3:0]
        vfnmadd231pd ymm10, ymm7, ymm2  #ymm10 = v.x - coil->R * cos_table[7:4]
        vfnmadd231pd ymm11, ymm7, ymm3  #ymm11 = v.z - coil->R * sin_table[7:4]

        vmulpd  ymm8, ymm4, ymm4
        vfmadd231pd ymm8, ymm5, ymm5
        vfmadd231pd ymm8, ymm6, ymm6    #ymm8 = Xr**2 + Yr**2 + Zr**2
        vmulpd  ymm12, ymm10, ymm10
        vfmadd231pd ymm12, ymm5, ymm5
        vfmadd231pd ymm12, ymm11, ymm11 #ymm12 = Xr**2 + Yr**2 + Zr**2

        vmulpd  ymm7, ymm0, ymm4
        vfmadd231pd ymm7, ymm1, ymm6    #Yc = Zj * Xr - Xj * Zr
        vmulpd  ymm0, ymm0, ymm5        #Xc = Zj * Yr
        vmulpd  ymm1, ymm1, ymm5        #Zc = -Xj * Yr
        vmulpd  ymm4, ymm2, ymm10
        vfmadd231pd ymm4, ymm3, ymm11   #Yc = Zj * Xr - Xj * Zr
        vmulpd  ymm2, ymm2, ymm5        #Xc = Zj * Yr
        vmulpd  ymm3, ymm3, ymm5        #Zc = -Xj * Yr

        #Reorder Vc, just for readability
        #0,7,1->0,1,2/2,4,3->3,4,5

        vmovapd ymm5, ymm3
        vmovapd ymm3, ymm2
        vmovapd ymm2, ymm1
        vmovapd ymm1, ymm7

        vcvtpd2ps xmm6, ymm8
        vrsqrtps xmm6, xmm6
        vcvtps2pd ymm6, xmm6            #Y0 = approx(1/sqrt(ymm8))
        vbroadcastsd ymm7, .onehalf[rip] #ymm7[3:0] = 1.5
        vmulpd  ymm8, ymm8, .neg_half[rip] #ymm8 = -0.5 * X

        vmulpd  ymm9, ymm8, ymm6        #ymm9 = -0.5 * X * Y0
        vfmadd213pd ymm9, ymm6, ymm7    #ymm9 = 1.5 - 0.5 * X * Y0 * Y0
        vmulpd  ymm6, ymm9, ymm6        #ymm6 = Y1 = Y0 * (1.5 - 0.5 * X * Y0 * Y0)
        vmulpd  ymm9, ymm8, ymm6        #ymm9 = -0.5 * X * Y1
        vfmadd213pd ymm9, ymm6, ymm7    #ymm9 = 1.5 - 0.5 * X * Y1 * Y1
        vmulpd  ymm6, ymm9, ymm6        #ymm6 = Y2 = Y1 * (1.5 - 0.5 * X * Y1 * Y1)
        vmulpd  ymm8, ymm6, ymm6
        vmulpd  ymm6, ymm6, ymm8        #ymm6 = 1.0/sqrt(X**2+Y**2+Z**2)**3

        vcvtpd2ps xmm13, ymm12
        vrsqrtps xmm13, xmm13
        vcvtps2pd ymm13, xmm13
        vmulpd  ymm12, ymm12, .neg_half[rip]

        vmulpd  ymm14, ymm13, ymm12
        vfmadd213pd ymm14, ymm13, ymm7
        vmulpd  ymm13, ymm14, ymm13
        vmulpd  ymm14, ymm13, ymm12
        vfmadd213pd ymm14, ymm13, ymm7
        vmulpd  ymm13, ymm14, ymm13
        vmulpd  ymm14, ymm13, ymm13
        vmulpd  ymm7, ymm13, ymm14      #ymm7 = 1.0/sqrt(X**2+Y**2+Z**2)**3

        vmulpd  ymm0, ymm0, ymm6
        vmulpd  ymm1, ymm1, ymm6
        vmulpd  ymm2, ymm2, ymm6
        vfmadd231pd ymm0, ymm3, ymm7
        vfmadd231pd ymm1, ymm4, ymm7
        vfmadd231pd ymm2, ymm5, ymm7

        vbroadcastsd ymm4, 32[rsp]
        vhaddpd ymm3, ymm0, ymm1
        vhaddpd ymm5, ymm2, ymm2
        vextractf128 xmm0, ymm3, 0x1
        vextractf128 xmm6, ymm5, 0x1
        vaddpd xmm6, xmm6, xmm5
        vaddpd xmm0, xmm0, xmm3
        vinsertf128 ymm0, ymm0, xmm6, 0x1

        vmulpd  ymm0, ymm0, ymm4
        vmovupd ymm1, [rdi]
        vaddpd  ymm0, ymm0, ymm1
        vmovupd [rdi], ymm0

        #vhaddpd ymm3, ymm0, ymm1
        #vhaddpd ymm4, ymm2, ymm2
        #vperm2f128 ymm5, ymm3, ymm4, 0x21
        #vblendpd ymm6, ymm3, ymm4, 0b1100
        #vaddpd  ymm7, ymm5, ymm6

        #vextractf128 xmm0, ymm7, 0x1
        #vaddpd  xmm7, xmm7, [rdi]
        #vmovapd [rdi], xmm7
        #movsd   xmm1, 16[rdi]
        #subsd   xmm1, xmm0
        #movsd   16[rdi], xmm1

        sub     edx, 1
        jne     biot_savart_loop1
        movsd   xmm0, 24[rsp]
        addsd   xmm0, 8[rsi]
        movsd   24[rsp], xmm0
        sub     ecx, 1
        jne     biot_savart_loop2
        addq    rsp, 40
        ret
        .size   biot_savart, .-biot_savart
        .section .rodata,"a",@progbits
        .align  8
.b_scale:
        .double 9.817477042468104e-09
.onehalf:
        .double 1.5
        .align  32
.neg_half:
        .double -0.5
        .double -0.5
        .double -0.5
        .double -0.5
.sin_table:
        .double 0.0
        .double 0.0980171403296
        .double 0.195090322016
        .double 0.290284677254
        .double 0.382683432365
        .double 0.471396736826
        .double 0.55557023302
        .double 0.634393284164
        .double 0.707106781187
        .double 0.773010453363
        .double 0.831469612303
        .double 0.881921264348
        .double 0.923879532511
        .double 0.956940335732
        .double 0.980785280403
        .double 0.995184726672
        .double 1.0
        .double 0.995184726672
        .double 0.980785280403
        .double 0.956940335732
        .double 0.923879532511
        .double 0.881921264348
        .double 0.831469612303
        .double 0.773010453363
        .double 0.707106781187
        .double 0.634393284164
        .double 0.55557023302
        .double 0.471396736826
        .double 0.382683432365
        .double 0.290284677254
        .double 0.195090322016
        .double 0.0980171403296
        .double 1.22464679915e-16
        .double -0.0980171403296
        .double -0.195090322016
        .double -0.290284677254
        .double -0.382683432365
        .double -0.471396736826
        .double -0.55557023302
        .double -0.634393284164
        .double -0.707106781187
        .double -0.773010453363
        .double -0.831469612303
        .double -0.881921264348
        .double -0.923879532511
        .double -0.956940335732
        .double -0.980785280403
        .double -0.995184726672
        .double -1.0
        .double -0.995184726672
        .double -0.980785280403
        .double -0.956940335732
        .double -0.923879532511
        .double -0.881921264348
        .double -0.831469612303
        .double -0.773010453363
        .double -0.707106781187
        .double -0.634393284164
        .double -0.55557023302
        .double -0.471396736826
        .double -0.382683432365
        .double -0.290284677254
        .double -0.195090322016
        .double -0.0980171403296
        .align  32
.cos_table:
        .double 1.0
        .double 0.995184726672
        .double 0.980785280403
        .double 0.956940335732
        .double 0.923879532511
        .double 0.881921264348
        .double 0.831469612303
        .double 0.773010453363
        .double 0.707106781187
        .double 0.634393284164
        .double 0.55557023302
        .double 0.471396736826
        .double 0.382683432365
        .double 0.290284677254
        .double 0.195090322016
        .double 0.0980171403296
        .double 6.12323399574e-17
        .double -0.0980171403296
        .double -0.195090322016
        .double -0.290284677254
        .double -0.382683432365
        .double -0.471396736826
        .double -0.55557023302
        .double -0.634393284164
        .double -0.707106781187
        .double -0.773010453363
        .double -0.831469612303
        .double -0.881921264348
        .double -0.923879532511
        .double -0.956940335732
        .double -0.980785280403
        .double -0.995184726672
        .double -1.0
        .double -0.995184726672
        .double -0.980785280403
        .double -0.956940335732
        .double -0.923879532511
        .double -0.881921264348
        .double -0.831469612303
        .double -0.773010453363
        .double -0.707106781187
        .double -0.634393284164
        .double -0.55557023302
        .double -0.471396736826
        .double -0.382683432365
        .double -0.290284677254
        .double -0.195090322016
        .double -0.0980171403296
        .double -1.83697019872e-16
        .double 0.0980171403296
        .double 0.195090322016
        .double 0.290284677254
        .double 0.382683432365
        .double 0.471396736826
        .double 0.55557023302
        .double 0.634393284164
        .double 0.707106781187
        .double 0.773010453363
        .double 0.831469612303
        .double 0.881921264348
        .double 0.923879532511
        .double 0.956940335732
        .double 0.980785280403
        .double 0.995184726672
