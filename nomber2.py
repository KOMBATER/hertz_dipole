import numpy as np
import matplotlib.pyplot as plt

f = 6.5e9
c = 3e8
lam = c / f
k = 2 * np.pi / lam
l = 0.005 * lam

theta = np.linspace(1e-6, np.pi - 1e-6, 2000)
theta_deg = theta * 180 / np.pi

E_theta = (np.cos(k * l * np.cos(theta)) - np.cos(k * l)) / np.sin(theta)

F = np.abs(E_theta)
F = F / np.max(F)

integral = np.trapz(F**2 * np.sin(theta), theta)
D_max = 4 * np.pi / (2 * np.pi * integral)
D_max_db = 10 * np.log10(D_max)

print("Максимальный КНД:")
print(f"D_max = {D_max:.4f}")
print(f"D_max = {D_max_db:.2f} дБ")

D_theta = D_max * F**2
D_theta_db = 10 * np.log10(D_theta)

file_name = "antenna_all_data.txt"

with open(file_name, "w", encoding="utf-8") as file:
    file.write("Таблица данных:\n")
    file.write("theta_rad    theta_deg    F(theta)    D(theta)    D(theta)_dB\n")

    for i in range(len(theta)):
        file.write(
            f"{theta[i]:.6e}  "
            f"{theta_deg[i]:.6f}  "
            f"{F[i]:.6e}  "
            f"{D_theta[i]:.6e}  "
            f"{D_theta_db[i]:.6f}\n"
        )

print(f"\nВсе данные сохранены в файл: {file_name}")

fig, axs = plt.subplots(2, 2, figsize=(12, 10))

axs[0, 0] = plt.subplot(2, 2, 1, polar=True)
axs[0, 0].plot(theta, D_theta_db)  # ДБ вместо разы
axs[0, 0].set_title("Полярная диаграмма КНД (дБ)")

axs[0, 1] = plt.subplot(2, 2, 2, polar=True)
axs[0, 1].plot(theta, D_theta)  # разы вместо ДБ
axs[0, 1].set_title("Полярная диаграмма КНД (разы)")

axs[1, 0].plot(theta_deg, D_theta_db)  # ДБ вместо разы
axs[1, 0].set_xlabel("θ, град")
axs[1, 0].set_ylabel("D(θ), дБ")
axs[1, 0].set_title("КНД в дБ")
axs[1, 0].grid()

axs[1, 1].plot(theta_deg, D_theta)  # разы вместо ДБ
axs[1, 1].set_xlabel("θ, град")
axs[1, 1].set_ylabel("D(θ)")
axs[1, 1].set_title("КНД в разах")
axs[1, 1].grid()

plt.tight_layout()
plt.show()