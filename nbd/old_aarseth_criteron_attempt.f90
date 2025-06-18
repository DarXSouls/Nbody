
      !ac3 = sqrt( axc3(i)*axc3(i) + ayc3(i)*ayc3(i) + azc3(i)*azc3(i) )
      !acx2_ext = axc2(i) + t_diff(i)*axc3(i)
      !acy2_ext = ayc2(i) + t_diff(i)*ayc3(i)
      !acz2_ext = azc2(i) + t_diff(i)*azc3(i)
      !ac2_ext = sqrt(acx2_ext*acx2_ext + acy2_ext*acy2_ext + acz2_ext*acz2_ext)
      !dt = min(dt_max, sqrt( eta_dt * (a1*ac2_ext + a1dot*a1dot) / (a1dot*ac3 + ac2_ext*ac2_ext) ))
