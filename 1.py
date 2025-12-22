
import numpy as np
from scipy.special import spherical_jn, spherical_yn
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from xml.dom import minidom


class EPRCalculator:
    def __init__(self, D: float, c=3e8, N: int = 100):
        self.c = c
        self.N = N
        self.r = D / 2

    def hankel_three_kind(self, n, k, r):
        return spherical_jn(n, k * r) + 1j * spherical_yn(n, k * r)

    def an(self, n, k, r):
        return spherical_jn(n, k * r) / self.hankel_three_kind(n, k, r)

    def bn(self, n, k, r):
        numerator = k * r * spherical_jn(n - 1, k * r) - n * spherical_jn(n, k * r)
        denominator = k * r * self.hankel_three_kind(n - 1, k, r) - n * self.hankel_three_kind(n, k, r)
        return numerator / denominator

    def sigma(self, lmbd, k, r):
        s = 0
        for n in range(1, self.N + 1):
            s += ((-1) ** n * (n + 0.5) * (self.bn(n, k, r) - self.an(n, k, r)))
        return (lmbd ** 2 / np.pi) * abs(s) ** 2

    def calculate(self, frequencies) -> list:
        sigmas = []
        for f in frequencies:
            lmbd = self.c / f
            k = 2 * np.pi / lmbd
            sigmas.append(self.sigma(lmbd, k, self.r))
        return sigmas


def save_to_xml(frequencies, sigmas, filename="result.xml"):
    root = ET.Element("data")

    for f, s in zip(frequencies, sigmas):
        row = ET.SubElement(root, "row")
        ET.SubElement(row, "freq").text = str(f)
        ET.SubElement(row, "lambda").text = str(3e8 / f)
        ET.SubElement(row, "rcs").text = str(s)

    rough_string = ET.tostring(root, encoding="utf-8")
    reparsed = minidom.parseString(rough_string)

    xml_declaration = '<?xml version="1.1" encoding="UTF-8"?>\n'
    pretty_xml = xml_declaration + reparsed.toprettyxml(indent="    ").split('\n', 1)[1]

    with open(filename, "w", encoding="utf-8") as f:
        f.write(pretty_xml)


if __name__ == '__main__':
    input_txt = 'task_rcs_02.txt'
    output_xml = 'result.xml'
    variant_numb = 1
    fmin = fmax = D = None
    with open(input_txt, encoding='utf-8') as file:
        lines = file.readlines()
    for line in lines:
        if line.strip() and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 4:
                variant = int(parts[0])
                if variant == variant_numb:
                    D = float(parts[1])
                    fmin = float(parts[2])
                    fmax = float(parts[3])
                    break

    frequencies = np.linspace(fmin, fmax, 200)

    calculator = EPRCalculator(D=D, N=50)
    sigmas = calculator.calculate(frequencies)

    save_to_xml(frequencies, sigmas, filename=output_xml)

    plt.plot(frequencies / 1e9, sigmas)
    plt.xlabel('Частота (ГГц)')
    plt.ylabel('ЭПР (м²)')
    plt.title(f'ЭПР сферы (D={D} м)')
    plt.grid(True)
    plt.savefig('rcs_plot.png')
    plt.show()
