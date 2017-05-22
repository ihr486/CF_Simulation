        .intel_syntax noprefix
        .text
        .p2align 4,,15
        .global biot_savart
biot_savart:
        sub     rsp, 40                 #Align stack to 16bytes
        vmovapd [rsp], xmm0
        vmovss  8[rsp], xmm1
        vzeroall
        vmovss  12[rsp], xmm0
        vmovss  xmm0, [rdi]
        vmulss  xmm0, xmm0, [rip + _b_scale]
        vmovss  36[rsp], xmm0           #32[rsp] = 1E-7*2*pi*coil->R/THETA_DIV
        mov     ecx, dword ptr 8[rdi]   #ecx = coil->N
        mov     eax, ecx
        sub     eax, 1
        vcvtsi2ss xmm0, xmm0, eax
        vmulss  xmm0, xmm0, [rip + _neg_half]
        vmulss  xmm0, xmm0, 4[rdi]
        vaddss  xmm0, xmm0, 4[rsp]
        vmovss  4[rsp], xmm0           #4[rsp] = v.y - 0.5 * coil->P * (coil->N - 1)
biot_savart_loop2:
        mov     edx, 4
        lea     r8, [rip + _cos_table]
        lea     r9, [rip + _sin_table]
biot_savart_loop1:
        vmovaps ymm0, [r8]
        vmovaps ymm2, 32[r8]
        vmovaps ymm1, [r9]
        vmovaps ymm3, 32[r9]
        add     r8, 64
        add     r9, 64

        vbroadcastss ymm4, [rsp]      #ymm4[3:0] = v.x
        vbroadcastss ymm5, 4[rsp]
        vbroadcastss ymm6, 8[rsp]      #ymm6[3:0] = v.z
        vbroadcastss ymm7, [rdi]        #ymm7[3:0] = coil->R
        vmovaps ymm10, ymm4
        vmovaps ymm11, ymm6

        vfnmadd231ps ymm4, ymm7, ymm0   #ymm4 = v.x - coil->R * cos_table[3:0]
        vfnmadd231ps ymm6, ymm7, ymm1   #ymm6 = v.z - coil->R * sin_table[3:0]
        vfnmadd231ps ymm10, ymm7, ymm2  #ymm10 = v.x - coil->R * cos_table[7:4]
        vfnmadd231ps ymm11, ymm7, ymm3  #ymm11 = v.z - coil->R * sin_table[7:4]

        vmulps  ymm8, ymm4, ymm4
        vfmadd231ps ymm8, ymm5, ymm5
        vfmadd231ps ymm8, ymm6, ymm6    #ymm8 = Xr**2 + Yr**2 + Zr**2
        vmulps  ymm12, ymm10, ymm10
        vfmadd231ps ymm12, ymm5, ymm5
        vfmadd231ps ymm12, ymm11, ymm11 #ymm12 = Xr**2 + Yr**2 + Zr**2

        vmulps  ymm7, ymm0, ymm4
        vfmadd231ps ymm7, ymm1, ymm6    #Yc = Zj * Xr - Xj * Zr
        vmulps  ymm0, ymm0, ymm5        #Xc = Zj * Yr
        vmulps  ymm1, ymm1, ymm5        #Zc = -Xj * Yr
        vmulps  ymm4, ymm2, ymm10
        vfmadd231ps ymm4, ymm3, ymm11   #Yc = Zj * Xr - Xj * Zr
        vmulps  ymm2, ymm2, ymm5        #Xc = Zj * Yr
        vmulps  ymm3, ymm3, ymm5        #Zc = -Xj * Yr

        #Reorder Vc, just for readability
        #0,7,1->0,1,2/2,4,3->3,4,5

        vmovaps ymm5, ymm3
        vmovaps ymm3, ymm2
        vmovaps ymm2, ymm1
        vmovaps ymm1, ymm7

        vrsqrtps ymm6, ymm8             #Y0 = approx(1/sqrt(ymm8))
        vbroadcastss ymm7, [rip + _three] #ymm7[3:0] = 3

        vmulps  ymm9, ymm8, ymm6        #ymm9 = X * Y0
        vfnmadd213ps ymm9, ymm6, ymm7   #ymm9 = 3 - X * Y0 * Y0
        vmulps  ymm6, ymm9, ymm6        #ymm6 = Y1 = Y0 * (3 - X * Y0 * Y0)
        vmulps  ymm8, ymm6, ymm6
        vmulps  ymm6, ymm6, ymm8        #ymm6 = 8.0/sqrt(X**2+Y**2+Z**2)**3

        vrsqrtps ymm13, ymm12

        vmulps  ymm14, ymm13, ymm12
        vfnmadd213ps ymm14, ymm13, ymm7
        vmulps  ymm13, ymm14, ymm13
        vmulps  ymm14, ymm13, ymm13
        vmulps  ymm7, ymm13, ymm14      #ymm7 = 8.0/sqrt(X**2+Y**2+Z**2)**3

        vdpps   ymm0, ymm0, ymm6, 0xFF
        vdpps   ymm1, ymm1, ymm6, 0xFF
        vdpps   ymm2, ymm2, ymm6, 0xFF
        vdpps   ymm3, ymm3, ymm7, 0xFF
        vdpps   ymm4, ymm4, ymm7, 0xFF
        vdpps   ymm5, ymm5, ymm7, 0xFF

        vaddps  ymm0, ymm0, ymm3
        vaddps  ymm1, ymm1, ymm4
        vaddps  ymm2, ymm2, ymm5

        vextractf128 xmm3, ymm0, 1
        vextractf128 xmm4, ymm1, 1
        vextractf128 xmm5, ymm2, 1
        vaddss  xmm0, xmm0, xmm3
        vaddss  xmm1, xmm1, xmm4
        vaddss  xmm2, xmm2, xmm5

        vunpcklps xmm3, xmm0, xmm1
        vshufps xmm0, xmm3, xmm2, 0x44

        vbroadcastss xmm3, 36[rsp]
        vfmadd231ps xmm15, xmm0, xmm3

        sub     edx, 1
        jne     biot_savart_loop1
        vmovss  xmm0, 4[rsp]
        vaddss  xmm0, xmm0, 4[rdi]
        vmovss  4[rsp], xmm0
        sub     ecx, 1
        jne     biot_savart_loop2
        vmovaps 16[rsp], xmm15
        vmovaps xmm0, xmm15
        vmovss  xmm1, 8[rsp]
        add     rsp, 40
        ret
        .data
        .align  4
_b_scale:
        .float  1.2271846303E-9
_neg_half:
        .float  -0.5
_three:
        .float  3
        .align  32
_sin_table:
        .float 0.0
        .float 0.0980171403296
        .float 0.195090322016
        .float 0.290284677254
        .float 0.382683432365
        .float 0.471396736826
        .float 0.55557023302
        .float 0.634393284164
        .float 0.707106781187
        .float 0.773010453363
        .float 0.831469612303
        .float 0.881921264348
        .float 0.923879532511
        .float 0.956940335732
        .float 0.980785280403
        .float 0.995184726672
        .float 1.0
        .float 0.995184726672
        .float 0.980785280403
        .float 0.956940335732
        .float 0.923879532511
        .float 0.881921264348
        .float 0.831469612303
        .float 0.773010453363
        .float 0.707106781187
        .float 0.634393284164
        .float 0.55557023302
        .float 0.471396736826
        .float 0.382683432365
        .float 0.290284677254
        .float 0.195090322016
        .float 0.0980171403296
        .float 1.22464679915e-16
        .float -0.0980171403296
        .float -0.195090322016
        .float -0.290284677254
        .float -0.382683432365
        .float -0.471396736826
        .float -0.55557023302
        .float -0.634393284164
        .float -0.707106781187
        .float -0.773010453363
        .float -0.831469612303
        .float -0.881921264348
        .float -0.923879532511
        .float -0.956940335732
        .float -0.980785280403
        .float -0.995184726672
        .float -1.0
        .float -0.995184726672
        .float -0.980785280403
        .float -0.956940335732
        .float -0.923879532511
        .float -0.881921264348
        .float -0.831469612303
        .float -0.773010453363
        .float -0.707106781187
        .float -0.634393284164
        .float -0.55557023302
        .float -0.471396736826
        .float -0.382683432365
        .float -0.290284677254
        .float -0.195090322016
        .float -0.0980171403296
_cos_table:
        .float 1.0
        .float 0.995184726672
        .float 0.980785280403
        .float 0.956940335732
        .float 0.923879532511
        .float 0.881921264348
        .float 0.831469612303
        .float 0.773010453363
        .float 0.707106781187
        .float 0.634393284164
        .float 0.55557023302
        .float 0.471396736826
        .float 0.382683432365
        .float 0.290284677254
        .float 0.195090322016
        .float 0.0980171403296
        .float 6.12323399574e-17
        .float -0.0980171403296
        .float -0.195090322016
        .float -0.290284677254
        .float -0.382683432365
        .float -0.471396736826
        .float -0.55557023302
        .float -0.634393284164
        .float -0.707106781187
        .float -0.773010453363
        .float -0.831469612303
        .float -0.881921264348
        .float -0.923879532511
        .float -0.956940335732
        .float -0.980785280403
        .float -0.995184726672
        .float -1.0
        .float -0.995184726672
        .float -0.980785280403
        .float -0.956940335732
        .float -0.923879532511
        .float -0.881921264348
        .float -0.831469612303
        .float -0.773010453363
        .float -0.707106781187
        .float -0.634393284164
        .float -0.55557023302
        .float -0.471396736826
        .float -0.382683432365
        .float -0.290284677254
        .float -0.195090322016
        .float -0.0980171403296
        .float -1.83697019872e-16
        .float 0.0980171403296
        .float 0.195090322016
        .float 0.290284677254
        .float 0.382683432365
        .float 0.471396736826
        .float 0.55557023302
        .float 0.634393284164
        .float 0.707106781187
        .float 0.773010453363
        .float 0.831469612303
        .float 0.881921264348
        .float 0.923879532511
        .float 0.956940335732
        .float 0.980785280403
        .float 0.995184726672
