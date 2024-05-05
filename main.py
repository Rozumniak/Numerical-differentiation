from sympy import symbols
import math

x = symbols('x')

y_list_values = [1.767, 2.011, 2.27, 2.545, 2.835, 3.142, 3.464]
x_values = [1.5, 2.1]
h = 0.1

x_list_values = [x_values[0]]
value = x_list_values[0]

for i in range(int((x_values[1] - x_values[0])/h)):
    value += h
    x_list_values.append(value)

def newton(x_value):
    dy1 = []
    dy2 = []
    dy3 = []
    dy4 = []
    dy5 = []
    dy6 = []
    dyy = [y_list_values, dy1, dy2, dy3, dy4, dy5, dy6]

    for i in range(0, len(y_list_values) - 1):
        for j in range(0, len(y_list_values) - i - 1):
            dyy[i + 1].append((dyy[i][j + 1] - dyy[i][j]))

    max_length = max(len(x_list_values), len(y_list_values), len(dy1), len(dy2), len(dy3), len(dy4))

    print("Значення кінцевих різниць")
    print("{:<8} | {:<10} | {:<10} | {:<10} | {:<10} | {:<10} | {:<10} | {:<10} |".format("x", "y", "dy1", "dy2", "dy3", "dy4", "dy5", "dy6", ))
    print("_____________________________________________________________________________________________________")

    for i in range(max_length):
        x_val = "{:<8.4f}".format(x_list_values[i]) if i < len(x_list_values) else "         "
        y_val = "{:<10.4f}".format(y_list_values[i]) if i < len(y_list_values) else "         "
        dy1_val = "{:<10.4f}".format(dy1[i]) if i < len(dy1) else "          "
        dy2_val = "{:<10.4f}".format(dy2[i]) if i < len(dy2) else "          "
        dy3_val = "{:<10.4f}".format(dy3[i]) if i < len(dy3) else "          "
        dy4_val = "{:<10.4f}".format(dy4[i]) if i < len(dy4) else "          "
        dy5_val = "{:<10.4f}".format(dy5[i]) if i < len(dy5) else "          "
        dy6_val = "{:<10.4f}".format(dy6[i]) if i < len(dy6) else "          "
        print("{} | {} | {} | {} | {} | {} | {} | {} |".format(x_val, y_val, dy1_val, dy2_val, dy3_val, dy4_val, dy5_val, dy6_val))

    t = symbols('t')
    sign = symbols('sign')

    del_y1, del_y2, del_y3, del_y4, del_y5, del_y6 = symbols('del_y1, del_y2, del_y3, del_y4, del_y5, del_y6')

    first_der = (del_y1 + ((2 * t + sign) / 2) * del_y2 +
                 ((3 * t ** 2 + (sign * 6 * t) + 2) / 6) * del_y3 -
                 ((4 * t ** 3 + (sign * 18 * t ** 2) + 22 * t + sign * 6) / 24) * del_y4 +
                 ((5 * t ** 4 + (sign * 40 * t ** 3) + 105 * t ** 2 + (sign * 100 * t) + 24) / 120) * del_y5 +
                 ((6 * t ** 5 + (sign * 75 * t ** 4) + 340 * t ** 3 + (
                             sign * 675 * t ** 2) + 548 * t + sign * 120) / 120) * del_y6) / h

    second_der = ((del_y2 + (t + sign) * del_y3
                   + ((6 * t ** 2 + (sign * 18 * t) + 11) / 12) * del_y4
                   + ((20 * t ** 3 + (sign * 120 * t ** 2) + 210 * t + sign * 100) / 120) * del_y5
                   + ((30 * t ** 4 + (sign * 300 * t ** 3) + 1020 * t ** 2 + sign * 1350 * t + 548) / 720) * del_y6)
                    / (h ** 2))

    analytic_1 = math.pi * x_value / 2
    analytic_2 = math.pi / 2

    print(f"Перша та друга похідні, обраховані аналітичним методом:\n "
          f"f'(x) = (Pi * x) / 2 ={analytic_1}\n"
          f"f''(x) = Pi / 2 ={analytic_2}\n")

    print(f"\n____Рішення за першою інтерполяційною формулою Ньютона____\n")

    t_value = (x_value - x_list_values[0]) / h
    newton_first_der = first_der.subs([(t, t_value),
                                       (del_y1, dy1[0]),
                                       (del_y2, dy2[0]),
                                       (del_y3, dy3[0]),
                                       (del_y4, dy4[0]),
                                       (del_y5, dy5[0]),
                                       (del_y6, dy6[0]),
                                       (sign, -1)]).evalf()
    newton_second_der = second_der.subs([(t, t_value),
                                       (del_y1, dy1[0]),
                                       (del_y2, dy2[0]),
                                       (del_y3, dy3[0]),
                                       (del_y4, dy4[0]),
                                       (del_y5, dy5[0]),
                                       (del_y6, dy6[0]),
                                       (sign, -1)]).evalf()

    print(f"Значення для f'(x) = {newton_first_der} за першою формулою Ньютона")
    delta_value = abs(analytic_1 - newton_first_der)
    print(f"∆ = |{analytic_1} - {newton_first_der}| = {delta_value}")
    print(f"δ = {delta_value} / {analytic_1} = {delta_value / analytic_1}")

    print(f"\nЗначення для f''(x) = {newton_second_der} за першою формулою Ньютона")
    delta_value = abs(analytic_2 - newton_second_der)
    print(f"∆ = |{analytic_2} - {newton_second_der}| = {delta_value}")
    print(f"δ = {delta_value} / {analytic_2} = {delta_value / analytic_2}")

    print(f"\n____Рішення за другою інтерполяційною формулою Ньютона____\n")

    t_value = (x_value - x_list_values[6]) / h
    newton_first_der = first_der.subs([(t, t_value),
                                       (del_y1, dy1[5]),
                                       (del_y2, dy2[4]),
                                       (del_y3, dy3[3]),
                                       (del_y4, dy4[2]),
                                       (del_y5, dy5[1]),
                                       (del_y6, dy6[0]),
                                       (sign, 1)]).evalf()
    newton_second_der = second_der.subs([(t, t_value),
                                         (del_y1, dy1[5]),
                                         (del_y2, dy2[4]),
                                         (del_y3, dy3[3]),
                                         (del_y4, dy4[2]),
                                         (del_y5, dy5[1]),
                                         (del_y6, dy6[0]),
                                         (sign, 1)]).evalf()

    print(f"Значення для f'(x) = {newton_first_der} за другою формулою Ньютона")
    delta_value = abs(analytic_1 - newton_first_der)
    print(f"∆ = |{analytic_1} - {newton_first_der}| = {delta_value}")
    print(f"δ = {delta_value} / {analytic_1} = {delta_value / analytic_1}")

    print(f"\nЗначення для f''(x) = {newton_second_der} за другою формулою Ньютона")
    delta_value = abs(analytic_2 - newton_second_der)
    print(f"∆ = |{analytic_2} - {newton_second_der}| = {delta_value}")
    print(f"δ = {delta_value} / {analytic_2} = {delta_value / analytic_2}")

def lagrange(x_value):
    print(f"\n____Рішення за формулою Лагранжа____\n")
    first_der_f1 = (y_list_values[3] - y_list_values[2]) / h
    first_der_f2 = (y_list_values[4] - y_list_values[2]) / (2 * h)
    second_der = (y_list_values[4] - 2 * y_list_values[3] + y_list_values[2])/ (h ** 2)

    analytic_1 = math.pi * x_value / 2
    analytic_2 = math.pi / 2

    print(f"Обраховане значеня за 1 формулою Лагранжа для f'(x) = {first_der_f1}")
    delta_value = abs(analytic_1 - first_der_f1)
    print(f"∆ = |{analytic_1} - {first_der_f1}| = {delta_value}")
    print(f"δ = {delta_value} / {analytic_1} = {delta_value / analytic_1}\n")

    print(f"Обраховане значеня за 2 формулою Лагранжа для f'(x) = {first_der_f2}")
    delta_value = abs(analytic_1 - first_der_f2)
    print(f"∆ = |{analytic_1} - {first_der_f2}| = {delta_value}")
    print(f"δ = {delta_value} / {analytic_1} = {delta_value / analytic_1}\n")

    print(f"Обраховане значеня за формулою Лагранжа для f''(x) = {second_der}")
    delta_value = abs(analytic_2 - second_der)
    print(f"∆ = |{analytic_2} - {second_der}| = {delta_value}")
    print(f"δ = {delta_value} / {analytic_2} = {delta_value / analytic_2}\n")



def main():
    print("Комп'ютерний практикум №7 \nВаріант №11 \nВиконав студент групи ПБ-21 \nРозумняк Руслан\n")
    x_value = 1.83
    newton(x_value)
    lagrange(x_value)

if __name__ == "__main__":
    main()
