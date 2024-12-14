import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk, Label, Button, StringVar, messagebox, Frame, Radiobutton, IntVar, Scale, HORIZONTAL, Entry
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def plot_sir_model(beta, gamma, N, S0, t_0, t_final, h, method, canvas, outbreak_label, a1_label, a2_label, a3_label):
    try:
        beta = float(beta)
        gamma = float(gamma)
        N = int(N)
        S0 = int(S0)
        t_0 = float(t_0)
        t_final = float(t_final)
        h = float(h)

        I0 = N - S0  # Initial infectious population
        R0 = 0       # Initial removed population

        # Calculate number of steps
        n = int(np.ceil((t_final - t_0) / h))

        # Initialize arrays
        t = np.linspace(t_0, t_final, n + 1)
        y = np.zeros((3, n + 1))  # Rows: S, I, R; Columns: time steps

        # Initial conditions
        y[:, 0] = [S0, I0, R0]

        # Define the SIR model
        def f(y):
            S, I, R = y
            dS = -beta * S * I
            dI = beta * S * I - gamma * I
            dR = gamma * I
            return np.array([dS, dI, dR])

        # Euler's method
        if method == 1:
            for i in range(1, n + 1):
                y[:, i] = y[:, i - 1] + h * f(y[:, i - 1])
            title = "Euler's Approximation"

        # Runge-Kutta 4th Order Method
        elif method == 2:
            for i in range(1, n + 1):
                k1 = h * f(y[:, i - 1])
                k2 = h * f(y[:, i - 1] + 0.5 * k1)
                k3 = h * f(y[:, i - 1] + 0.5 * k2)
                k4 = h * f(y[:, i - 1] + k3)
                y[:, i] = y[:, i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            title = "Runge-Kutta Approximation"

        # Calculate A1, A2, A3
        dS_dt = -beta * y[0, :] * y[1, :]
        A1 = t[np.argmin(dS_dt)]
        A2_idx = np.argmax(y[1, :])
        A2 = (t[A2_idx], (y[1, A2_idx] / N) * 100)

        # Determine outbreak status and calculate A3 accordingly
        peak_reached = False
        outbreak_ended_day = None
        for i in range(1, len(y[1, :])):
            if y[1, i] < y[1, i - 1]:
                peak_reached = True
            if peak_reached and y[1, i] < 0.001 * N:
                outbreak_status = f"Outbreak Ended at t = {t[i]:.2f}"
                outbreak_ended_day = i
                break
        else:
            outbreak_status = "Outbreak Ongoing"

        if outbreak_ended_day is not None:
            A3 = (y[2, outbreak_ended_day] / N) * 100
        else:
            A3 = (y[2, -1] / N) * 100

        # Update outbreak status label
        outbreak_label.config(text=f"{outbreak_status}")

        # Update A1, A2, A3 labels
        a1_label.config(text=f"A1: Disease spreading most rapidly at t = {A1:.2f}")
        a2_label.config(text=f"A2: Peak on day {A2[0]:.2f}, {A2[1]:.2f}% infectious")
        a3_label.config(text=f"A3: {A3:.2f}% of population infected")

        # Create the plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(t, y[0, :], label='S(t)', linewidth=2)
        ax.plot(t, y[1, :], label='I(t)', linewidth=2)
        ax.plot(t, y[2, :], label='R(t)', linewidth=2)

        ax.set_title(title, fontsize=14)
        ax.set_xlabel("Time (t)", fontsize=12)
        ax.set_ylabel("Population", fontsize=12)
        ax.legend(fontsize=12)
        ax.grid()

        # Clear the canvas and draw the new plot
        for widget in canvas.winfo_children():
            widget.destroy()

        canvas_plot = FigureCanvasTkAgg(fig, master=canvas)
        canvas_plot.draw()
        canvas_plot.get_tk_widget().pack(fill='both', expand=True)

        plt.close(fig)

    except ValueError as e:
        messagebox.showerror("Invalid Input", "Please ensure all inputs are valid numbers.")

# Create the GUI
def main():
    def update_plot(*args):
        plot_sir_model(
            sliders["Beta"].get(),
            sliders["Gamma"].get(),
            n_var.get(),
            s0_var.get(),
            t0_var.get(),
            tfinal_var.get(),
            sliders["h"].get(),
            method.get(),
            canvas,
            outbreak_label,
            a1_label,
            a2_label,
            a3_label
        )

    root = Tk()
    root.title("SIR Model GUI")
    root.geometry("1200x700")

    # Input labels and sliders
    inputs = {
        "Beta (Infection rate)": 5e-5,
        "Gamma (Removal rate)": 0.142857,
        "N (Total population)": 50000,
        "S0 (Initial susceptible population)": 49990,
        "t_0 (Initial time)": 0,
        "t_final (Final time)": 53,
        "h (Time step)": 0.2
    }

    frame_inputs = Frame(root)
    frame_inputs.pack(side="right", padx=10, pady=10)

    sliders = {}

    Label(frame_inputs, text="Beta (Infection rate):").grid(row=0, column=0, sticky="e", pady=5)
    sliders["Beta"] = Scale(frame_inputs, from_=0.00001, to=0.001, resolution=0.000001, orient=HORIZONTAL, length=200, command=update_plot)
    sliders["Beta"].set(inputs["Beta (Infection rate)"])
    sliders["Beta"].grid(row=0, column=1, pady=5)

    Label(frame_inputs, text="Gamma (Removal rate):").grid(row=1, column=0, sticky="e", pady=5)
    sliders["Gamma"] = Scale(frame_inputs, from_=0.001, to=1.0, resolution=0.001, orient=HORIZONTAL, length=200, command=update_plot)
    sliders["Gamma"].set(inputs["Gamma (Removal rate)"])
    sliders["Gamma"].grid(row=1, column=1, pady=5)

    Label(frame_inputs, text="h (Time step):").grid(row=2, column=0, sticky="e", pady=5)
    sliders["h"] = Scale(frame_inputs, from_=0.1, to=1.0, resolution=0.1, orient=HORIZONTAL, length=200, command=update_plot)
    sliders["h"].set(inputs["h (Time step)"])
    sliders["h"].grid(row=2, column=1, pady=5)

    Label(frame_inputs, text="N (Total population):").grid(row=3, column=0, sticky="e", pady=5)
    n_var = StringVar(value=str(inputs["N (Total population)"]))
    Entry(frame_inputs, textvariable=n_var, width=20).grid(row=3, column=1, pady=5)

    Label(frame_inputs, text="S0 (Initial susceptible population):").grid(row=4, column=0, sticky="e", pady=5)
    s0_var = StringVar(value=str(inputs["S0 (Initial susceptible population)"]))
    Entry(frame_inputs, textvariable=s0_var, width=20).grid(row=4, column=1, pady=5)

    Label(frame_inputs, text="t_0 (Initial time):").grid(row=5, column=0, sticky="e", pady=5)
    t0_var = StringVar(value=str(inputs["t_0 (Initial time)"]))
    Entry(frame_inputs, textvariable=t0_var, width=20).grid(row=5, column=1, pady=5)

    Label(frame_inputs, text="t_final (Final time):").grid(row=6, column=0, sticky="e", pady=5)
    tfinal_var = StringVar(value=str(inputs["t_final (Final time)"]))
    Entry(frame_inputs, textvariable=tfinal_var, width=20).grid(row=6, column=1, pady=5)

    # Frame for A1, A2, A3 values
    frame_metrics = Frame(root)
    frame_metrics.pack(side="left", padx=10, pady=10)

    outbreak_label = Label(frame_metrics, text="Outbreak Status: N/A", font=("Arial", 14), fg="blue")
    outbreak_label.pack(anchor="w", pady=5)

    a1_label = Label(frame_metrics, text="A1: Disease spreading most rapidly at t = N/A", font=("Arial", 12))
    a1_label.pack(anchor="w", pady=5)

    a2_label = Label(frame_metrics, text="A2: Peak on day N/A, N/A% infectious", font=("Arial", 12))
    a2_label.pack(anchor="w", pady=5)

    a3_label = Label(frame_metrics, text="A3: N/A% of population infected", font=("Arial", 12))
    a3_label.pack(anchor="w", pady=5)

    # Canvas for plot
    canvas = Frame(root, width=700, height=600, bg="white")
    canvas.pack(side="left", padx=10, pady=10, fill="both", expand=True)

    # Radio buttons for method selection
    method = IntVar(value=1)
    Radiobutton(frame_inputs, text="Euler", variable=method, value=1).grid(row=7, column=0, pady=5, sticky="w")
    Radiobutton(frame_inputs, text="Runge-Kutta", variable=method, value=2).grid(row=7, column=1, pady=5, sticky="w")

    # Submit button
    Button(frame_inputs, text="Submit", command=lambda: plot_sir_model(
        sliders["Beta"].get(),
        sliders["Gamma"].get(),
        n_var.get(),
        s0_var.get(),
        t0_var.get(),
        tfinal_var.get(),
        sliders["h"].get(),
        method.get(),
        canvas,
        outbreak_label,
        a1_label,
        a2_label,
        a3_label
    )).grid(row=8, column=0, columnspan=1, pady=20)

    # Close button
    Button(frame_inputs, text="Close", command=root.quit).grid(row=8, column=1, columnspan=1, pady=20)

    # Run the GUI
    root.mainloop()

if __name__ == "__main__":
    main()
