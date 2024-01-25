def calc_uinf_ub(delta, mode, dir):
    if dir == "left":
        norm = -1.0
        uinf = norm * delta
        u_L = uinf - 1e-14
    elif dir == "right":
        norm = 1.0
        uinf = norm * delta
        u_L = uinf
    else:
        raise ValueError(f"invalid dir {dir}")

    if mode == "code":
        return uinf, u_L + (uinf - u_L * uinf * norm) * norm
    elif mode == "email":
        return uinf, u_L + (uinf - u_L * norm) * norm
    elif mode == "op":
        return uinf, norm * (u_L + (uinf - u_L * norm) * norm)
    elif mode == "uinf":
        return uinf, uinf
    else:
        raise ValueError(f"invalid mode {mode}")


def main(delta=1.0):
    print(f"delta={delta}")
    for dir in ["left", "right"]:
        s = ""
        for mode in ["code", "email", "op"]:
            uinf, val = calc_uinf_ub(delta, mode, dir)
            if not s:
                s = f" {dir}:  uinf={uinf}"
            s += f"  {mode}={val}"
        print(s)
    print()


main(0.7)
main(1.0)
