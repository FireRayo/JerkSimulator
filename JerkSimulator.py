#JerkSimulator Version 1.0

"""
Reference:
https://hidraulica.fluidas.ro/wp-content/uploads/07-25.pdf

"""

import sys
import math
import time
import csv
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from pathlib import Path
import sys, tkinter as tk


# -------------------- Iconos --------------------

def resource_path(rel):
    base = getattr(sys, "_MEIPASS", Path(__file__).resolve().parent)
    return str(Path(base, rel))

class SuiteApp(tk.Tk):
    def __init__(self):
        super().__init__()
        try:
            img = tk.PhotoImage(file=resource_path("assets/app.png"))
            self.iconphoto(True, img)
            self._icon_keep = img
        except Exception:
            pass

# -------------------- Núcleo de perfil JerkSimulator --------------------

def _durations_for_peak_velocity(Vp, a_max, j_max):
    if Vp <= 0 or a_max <= 0 or j_max <= 0:
        return 0.0, 0.0, True
    V_thr = (a_max**2) / j_max
    if Vp <= V_thr:
        t_j = math.sqrt(max(Vp, 0.0) / j_max)
        t_a_hold = 0.0
        triangular = True
    else:
        t_j = a_max / j_max
        V_from_jerk = (a_max**2) / j_max
        t_a_hold = max(0.0, (Vp - V_from_jerk) / a_max)
        triangular = False
    return t_j, t_a_hold, triangular


def _integrate_segments(a_sign, t_j, t_hold, j_max, dt, t_list, j_list, a_list, v_list, x_list):
    t = t_list[-1]; a = a_list[-1]; v = v_list[-1]; x = x_list[-1]

    def step(j_cmd, duration):
        nonlocal t, a, v, x
        if duration <= 0.0:
            return
        steps = max(1, int(round(duration / dt)))
        h = duration / steps
        for _ in range(steps):
            a = a + j_cmd * h
            v = v + a * h
            x = x + v * h
            t = t + h
            t_list.append(t); j_list.append(j_cmd); a_list.append(a); v_list.append(v); x_list.append(x)

    step(a_sign * j_max, t_j)
    step(0.0, t_hold)
    step(-a_sign * j_max, t_j)


def _build_0_to_Vp_to_0(Vp, a_max, d_max, j_max, dt):
    t_list = [0.0]; j_list = [0.0]; a_list = [0.0]; v_list = [0.0]; x_list = [0.0]
    if (Vp <= 0.0) or (a_max <= 0.0) or (d_max <= 0.0) or (j_max <= 0.0):
        return t_list, j_list, a_list, v_list, x_list
    t_j_acc, t_ah, _ = _durations_for_peak_velocity(Vp, a_max, j_max)
    _integrate_segments(+1.0, t_j_acc, t_ah, j_max, dt, t_list, j_list, a_list, v_list, x_list)
    t_j_dec, t_dh, _ = _durations_for_peak_velocity(Vp, d_max, j_max)
    _integrate_segments(-1.0, t_j_dec, t_dh, j_max, dt, t_list, j_list, a_list, v_list, x_list)
    return t_list, j_list, a_list, v_list, x_list


def distance_for_peak_velocity(Vp, a_max, d_max, j_max, dt):
    t, j, a, v, x = _build_0_to_Vp_to_0(Vp, a_max, d_max, j_max, dt)
    return x[-1] if len(x) else 0.0


def build_move_profile(distance, V_max, a_max, d_max, j_max, dt):
    t_list = [0.0]; j_list = [0.0]; a_list = [0.0]; v_list = [0.0]; x_list = [0.0]
    if distance <= 0.0 or V_max <= 0.0 or a_max <= 0.0 or d_max <= 0.0 or j_max <= 0.0:
        return {"t": np.array(t_list), "v": np.array(v_list), "a": np.array(a_list),
                "j": np.array(j_list), "x": np.array(x_list)}

    s_min = distance_for_peak_velocity(V_max, a_max, d_max, j_max, dt=dt)

    if distance >= s_min + 1e-9:
        t_j_acc, t_ah, _ = _durations_for_peak_velocity(V_max, a_max, j_max)
        t_j_dec, t_dh, _ = _durations_for_peak_velocity(V_max, d_max, j_max)
        _integrate_segments(+1.0, t_j_acc, t_ah, j_max, dt, t_list, j_list, a_list, v_list, x_list)
        s_cruise = max(0.0, distance - s_min)
        if V_max > 0 and s_cruise > 0.0:
            t_cruise = s_cruise / V_max
            steps = max(1, int(round(t_cruise / dt)))
            h = t_cruise / steps
            t = t_list[-1]; x = x_list[-1]; v = V_max
            for _ in range(steps):
                x += v * h; t += h
                t_list.append(t); j_list.append(0.0); a_list.append(0.0); v_list.append(v); x_list.append(x)
        _integrate_segments(-1.0, t_j_dec, t_dh, j_max, dt, t_list, j_list, a_list, v_list, x_list)
        return {"t": np.array(t_list), "v": np.array(v_list), "a": np.array(a_list),
                "j": np.array(j_list), "x": np.array(x_list)}

    lo, hi = 0.0, V_max
    target = distance
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        s_mid = distance_for_peak_velocity(mid, a_max, d_max, j_max, dt=dt)
        if s_mid < target: lo = mid
        else: hi = mid
        if abs(s_mid - target) <= max(1e-3 * target, 0.1): break
    V_peak = max(0.0, min(V_max, 0.5 * (lo + hi)))
    t, j, a, v, x = _build_0_to_Vp_to_0(V_peak, a_max, d_max, j_max, dt)
    if len(x):
        dx = distance - x[-1]
        if abs(dx) > 1e-9 and V_peak > 0.0:
            t_patch = dx / max(V_peak, 1e-9)
            t0 = t[-1]; x0 = x[-1]; v0 = v[-1]
            steps = max(1, int(round(abs(t_patch) / dt)))
            h = (t_patch / steps) if steps > 0 else 0.0
            for _ in range(max(0, steps)):
                t0 += h; x0 += v0 * h
                t.append(t0); j.append(0.0); a.append(0.0); v.append(v0); x.append(x0)
    return {"t": np.array(t), "v": np.array(v), "a": np.array(a),
            "j": np.array(j), "x": np.array(x)}


# ------------------------------ GUI con MAX + % ------------------------------

class SuiteApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("JerkSimulator V1.0")
        self.geometry("1250x900"); self.minsize(1050, 760)

        # Estado animación
        self._pending_after = None
        self._timer = None
        self._playing = False
        self.play_time = 0.0
        self._last_tick = None
        self.speed = 1.0
        self.loop = tk.BooleanVar(value=False)
        self.fps = tk.IntVar(value=60)
        self._frame_ms = 16

        # Máximos editables
        self.vmax_max = tk.DoubleVar(value=400.0)
        self.amax_max = tk.DoubleVar(value=500.0)
        self.dmax_max = tk.DoubleVar(value=500.0)
        self.jmax_max = tk.DoubleVar(value=1000.0)
        self.dist_max = tk.DoubleVar(value=2000.0)

        # Porcentajes (0-100)
        self.vmax_pct = tk.DoubleVar(value=50.0)
        self.amax_pct = tk.DoubleVar(value=50.0)
        self.dmax_pct = tk.DoubleVar(value=50.0)
        self.jmax_pct = tk.DoubleVar(value=50.0)
        self.dist_pct = tk.DoubleVar(value=50.0)

        # Δt de integración
        self.dtms = tk.DoubleVar(value=1.0)

        # Datos resultantes
        self.data = {"t": np.array([0.0]), "v": np.array([0.0]), "a": np.array([0.0]),
                     "j": np.array([0.0]), "x": np.array([0.0])}
        self._xmax_canvas = 1.0

        # Readouts para mostrar valor efectivo (texto a la derecha de cada slider)
        self._read_v = tk.StringVar(); self._read_a = tk.StringVar()
        self._read_d = tk.StringVar(); self._read_j = tk.StringVar(); self._read_x = tk.StringVar()

        self._controls = []
        self._build_ui()
        self._wire_traces()
        self._refresh_readouts()
        self.update_all()

    # ---------- Helpers UI ----------
    def _add_control(self, w):
        self._controls.append(w); return w

    def _set_controls_enabled(self, enabled: bool):
        for w in self._controls:
            try:
                if hasattr(w, "state"):
                    w.state(["!disabled"] if enabled else ["disabled"])
                else:
                    w.configure(state=("normal" if enabled else "disabled"))
            except Exception:
                pass

    def _wire_traces(self):
        # Cuando cambian máximos, refrescamos readouts y (si no está jugando) recalculamos
        for var in (self.vmax_max, self.amax_max, self.dmax_max, self.jmax_max, self.dist_max):
            var.trace_add("write", lambda *_: (self._refresh_readouts(), self.schedule_update()))
        # Cuando cambian %
        for var in (self.vmax_pct, self.amax_pct, self.dmax_pct, self.jmax_pct, self.dist_pct):
            var.trace_add("write", lambda *_: (self._refresh_readouts(), self.schedule_update()))

    def _effective_params(self):
        vmax = max(0.0, float(self.vmax_max.get()) * float(self.vmax_pct.get()) / 100.0)
        amax = max(0.0, float(self.amax_max.get()) * float(self.amax_pct.get()) / 100.0)
        dmax = max(0.0, float(self.dmax_max.get()) * float(self.dmax_pct.get()) / 100.0)
        jmax = max(0.0, float(self.jmax_max.get()) * float(self.jmax_pct.get()) / 100.0)
        dist = max(0.0, float(self.dist_max.get()) * float(self.dist_pct.get()) / 100.0)
        return vmax, amax, dmax, jmax, dist

    def _refresh_readouts(self):
        vmax, amax, dmax, jmax, dist = self._effective_params()
        self._read_v.set(f"{vmax:0.0f} mm/s  ({self.vmax_pct.get():0.0f}%)")
        self._read_a.set(f"{amax:0.0f} mm/s² ({self.amax_pct.get():0.0f}%)")
        self._read_d.set(f"{dmax:0.0f} mm/s² ({self.dmax_pct.get():0.0f}%)")
        self._read_j.set(f"{jmax:0.0f} mm/s³ ({self.jmax_pct.get():0.0f}%)")
        self._read_x.set(f"{dist:0.0f} mm    ({self.dist_pct.get():0.0f}%)")

    # ---------- UI Layout ----------
    def _build_ui(self):
        root = ttk.Frame(self); root.pack(side=tk.TOP, fill=tk.X, padx=10, pady=8)

        
        root.columnconfigure(2, weight=1)
# Encabezados de columnas
        hdr = ["Magnitude", "Maximum", "(0–100) %", "Current Value"]
        for i, h in enumerate(hdr):
            sticky = "ew" if i == 2 else "w"
            anchor = "center" if i == 2 else "w"
            ttk.Label(root, text=h, font=("TkDefaultFont", 12, "bold"), anchor=anchor)\
                .grid(row=0, column=i, sticky=sticky, padx=6, pady=(0,4))
# Filas con Entry de max + slider % + readout
        def row(idx, name, max_var, pct_var, read_var, unit):
            ttk.Label(root, text=f"{name} ({unit})").grid(row=idx, column=0, sticky="w", padx=6, pady=3)
            self._add_control(ttk.Entry(root, textvariable=max_var, width=10)).grid(row=idx, column=1, sticky="w", padx=6)
            s = self._add_control(ttk.Scale(root, from_=0, to=100, orient=tk.HORIZONTAL, variable=pct_var))
            s.grid(row=idx, column=2, sticky="ew", padx=6)
            root.columnconfigure(2, weight=1)
            ttk.Label(root, textvariable=read_var, width=20).grid(row=idx, column=3, sticky="e", padx=6)

        row(1, "Speed", self.vmax_max, self.vmax_pct, self._read_v, "mm/s")
        row(2, "Acceleration", self.amax_max, self.amax_pct, self._read_a, "mm/s²")
        row(3, "Deceleration", self.dmax_max, self.dmax_pct, self._read_d, "mm/s²")
        row(4, "Jerk", self.jmax_max, self.jmax_pct, self._read_j, "mm/s³")
        row(5, "Distance", self.dist_max, self.dist_pct, self._read_x, "mm")

        # Δt, Botones, FPS y Loop
        ctl = ttk.Frame(self); ctl.pack(side=tk.TOP, fill=tk.X, padx=10, pady=(6,8))
        ttk.Label(ctl, text="Δt (ms)").pack(side=tk.LEFT, padx=(0,6))
        self._add_control(ttk.Entry(ctl, textvariable=self.dtms, width=7)).pack(side=tk.LEFT, padx=(0,10))
        self._add_control(ttk.Button(ctl, text="Calculate", command=self.update_all)).pack(side=tk.LEFT, padx=4)
        self._add_control(ttk.Button(ctl, text="Export CSV", command=self.export_csv)).pack(side=tk.LEFT, padx=4)
        ttk.Label(ctl, text="FPS").pack(side=tk.LEFT, padx=(16,4))
        self.fps_box = ttk.Combobox(ctl, values=["30","60","90","120","144"], width=4, state="readonly")
        self.fps_box.current(1); self.fps_box.pack(side=tk.LEFT, padx=4)
        self.fps_box.bind("<<ComboboxSelected>>", self._on_fps_change)
        self._add_control(ttk.Checkbutton(ctl, text="Loop", variable=self.loop)).pack(side=tk.LEFT, padx=10)

        # Notebook con gráficas y animación
        notebook = ttk.Notebook(self); notebook.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=(0,10))

        # Gráficas
        tab_plot = ttk.Frame(notebook); notebook.add(tab_plot, text="Graphics")
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax_v = self.fig.add_subplot(411)
        self.ax_a = self.fig.add_subplot(412, sharex=self.ax_v)
        self.ax_j = self.fig.add_subplot(413, sharex=self.ax_v)
        self.ax_x = self.fig.add_subplot(414, sharex=self.ax_v)
        for ax, lab in zip((self.ax_v, self.ax_a, self.ax_j, self.ax_x),
                           ("v [mm/s]", "a [mm/s²]", "j [mm/s³]", "x [mm]")):
            ax.set_ylabel(lab); ax.grid(True, linestyle="--", alpha=0.4)
        self.ax_x.set_xlabel("t [s]")
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=tab_plot)
        self.canvas_plot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Animación
        tab_anim = ttk.Frame(notebook); notebook.add(tab_anim, text="Animation")
        top = ttk.Frame(tab_anim); top.pack(side=tk.TOP, fill=tk.X, padx=6, pady=6)
        self.play_btn = ttk.Button(top, text="Play", command=self.toggle_play); self.play_btn.pack(side=tk.LEFT, padx=4)
        ttk.Button(top, text="Reset", command=self.reset_anim).pack(side=tk.LEFT, padx=4)
        ttk.Label(top, text="Playback speed").pack(side=tk.LEFT, padx=(12,4))
        self.speed_box = ttk.Combobox(top, values=["0.25x","0.5x","1x","2x","4x"], width=6, state="readonly")
        self.speed_box.current(2); self.speed_box.pack(side=tk.LEFT, padx=4); self.speed_box.bind("<<ComboboxSelected>>", self._on_speed)

        self.canvas_w = 1100; self.canvas_h = 280
        self.canvas = tk.Canvas(tab_anim, width=self.canvas_w, height=self.canvas_h, bg="#101010", highlightthickness=0)
        self.canvas.pack(side=tk.TOP, fill=tk.X, expand=False, padx=10, pady=4)
        self.margin = 40
        y = self.canvas_h // 2
        self.canvas.create_line(self.margin, y, self.canvas_w - self.margin, y, width=4, tags="track")
        self.canvas.create_line(self.margin, y-12, self.margin, y+12, width=2, tags="tick_l")
        self.canvas.create_line(self.canvas_w - self.margin, y-12, self.canvas_w - self.margin, y+12, width=2, tags="tick_r")
        self.marker_r = 10
        self.marker_id = self.canvas.create_oval(0,0,0,0, outline="", fill="#e0e0e0")
        self._place_marker(0.0)

        info = ttk.Frame(tab_anim); info.pack(side=tk.TOP, fill=tk.X, padx=10, pady=6)
        self.lbl_t = ttk.Label(info, text="t: 0.000 s"); self.lbl_t.pack(side=tk.LEFT, padx=8)
        self.lbl_x = ttk.Label(info, text="x: 0.000 mm"); self.lbl_x.pack(side=tk.LEFT, padx=8)
        self.lbl_v = ttk.Label(info, text="v: 0.000 mm/s"); self.lbl_v.pack(side=tk.LEFT, padx=8)
        self.lbl_a = ttk.Label(info, text="a: 0.000 mm/s²"); self.lbl_a.pack(side=tk.LEFT, padx=8)
        self.lbl_j = ttk.Label(info, text="j: 0.000 mm/s³"); self.lbl_j.pack(side=tk.LEFT, padx=8)

        self._on_fps_change()  # init

    # -------------- Lógica --------------
    def schedule_update(self):
        if getattr(self, "_playing", False):
            return
        if self._pending_after is not None:
            try: self.after_cancel(self._pending_after)
            except Exception: pass
        self._pending_after = self.after(120, self.update_all)

    def update_all(self):
        self._playing = False; self._last_tick = None
        if self._timer is not None:
            try: self.after_cancel(self._timer)
            except Exception: pass
            self._timer = None
        if hasattr(self, "play_btn"):
            self.play_btn.config(text="Play")
        self._pending_after = None

        try:
            vmax, amax, dmax, jmax, dist = self._effective_params()
            dt = max(0.0001, float(self.dtms.get())/1000.0)
        except ValueError:
            messagebox.showerror("Error", "Invalid values ​​in controls."); return

        self.data = build_move_profile(dist, vmax, amax, dmax, jmax, dt=dt)

        # Precalcular x_max para canvas
        x = self.data["x"]
        self._xmax_canvas = float(np.max(x)) if len(x) else 1.0
        if self._xmax_canvas <= 1e-12: self._xmax_canvas = 1.0

        self._update_plots()
        self.play_time = 0.0
        self._update_anim_at_time(0.0)

    def _update_plots(self):
        t = self.data["t"]; v = self.data["v"]; a = self.data["a"]; j = self.data["j"]; x = self.data["x"]
        for ax in (self.ax_v, self.ax_a, self.ax_j, self.ax_x): ax.cla()
        self.ax_v.plot(t, v, linewidth=1.5); self.ax_v.set_ylabel("Speed\n[mm/s]")
        self.ax_a.plot(t, a, linewidth=1.2); self.ax_a.set_ylabel("Acceleration\nDeceleration\n[mm/s²]")
        self.ax_j.plot(t, j, linewidth=1.0); self.ax_j.set_ylabel("Jerk\n[mm/s³]")
        self.ax_x.plot(t, x, linewidth=1.5); self.ax_x.set_ylabel("Distance\n[mm]"); self.ax_x.set_xlabel("Time [s]")
        if len(t) > 1: self.ax_v.set_xlim(t[0], t[-1])
        for ax in (self.ax_v, self.ax_a, self.ax_j, self.ax_x): ax.grid(True, linestyle="--", alpha=0.4)
        self.fig.tight_layout(); self.canvas_plot.draw_idle()

    # -------- Animación --------
    def toggle_play(self):
        if self._pending_after is not None:
            try: self.after_cancel(self._pending_after)
            except Exception: pass
            self._pending_after = None

        if not self._playing:
            if len(self.data["t"]) <= 1: self.update_all()
            self.play_time = 0.0
            self._last_tick = time.perf_counter()
            self._playing = True; self.play_btn.config(text="Pause")
            self._set_controls_enabled(False)
            self._schedule_tick()
        else:
            self._playing = False; self.play_btn.config(text="Play")
            self._set_controls_enabled(True)
            if self._timer is not None:
                try: self.after_cancel(self._timer)
                except Exception: pass
                self._timer = None

    def reset_anim(self):
        self._playing = False; self.play_btn.config(text="Play")
        self._set_controls_enabled(True); self.play_time = 0.0; self._last_tick = None
        if self._timer is not None:
            try: self.after_cancel(self._timer)
            except Exception: pass
            self._timer = None
        self._update_anim_at_time(0.0)

    def _on_speed(self, *_):
        val = self.speed_box.get().replace("x","")
        try: self.speed = float(val)
        except: self.speed = 1.0

    def _on_fps_change(self, *_):
        try: fps = int(self.fps_box.get())
        except Exception: fps = 60
        self.fps.set(max(10, min(240, fps)))
        self._frame_ms = max(1, int(round(1000.0 / self.fps.get())))

    def _schedule_tick(self):
        if self._timer is not None:
            try: self.after_cancel(self._timer)
            except Exception: pass
        self._timer = self.after(self._frame_ms, self._tick)

    def _tick(self):
        if not self._playing:
            return
        now = time.perf_counter()
        dt = 0.016 if self._last_tick is None else (now - self._last_tick)
        self._last_tick = now
        self.play_time += dt * self.speed

        t = self.data["t"]
        if len(t) >= 2:
            total = t[-1]
            if total <= 0:
                self.play_time = 0.0
            else:
                if self.loop.get():
                    if self.play_time >= total:
                        self.play_time = self.play_time % total
                        self._update_anim_at_time(0.0)
                else:
                    if self.play_time >= total:
                        self.play_time = total
                        self._playing = False
                        self.play_btn.config(text="Play")
                        self._set_controls_enabled(True)
                        self._update_anim_at_time(self.play_time)
                        return

        self._update_anim_at_time(self.play_time)
        self._schedule_tick()

    def _update_anim_at_time(self, t_now):
        t = self.data["t"]; x = self.data["x"]; v = self.data["v"]; a = self.data["a"]; j = self.data["j"]
        if len(t) == 0:
            self._place_marker(0.0); self._update_labels(0.0, 0.0, 0.0, 0.0, 0.0); return

        if t_now <= t[0]:
            xi, vi, ai, ji = x[0], v[0], a[0], j[0]
            self._place_marker(xi); self._update_labels(t[0], xi, vi, ai, ji); return
        if t_now >= t[-1]:
            xi, vi, ai, ji = x[-1], v[-1], a[-1], j[-1]
            self._place_marker(xi); self._update_labels(t[-1], xi, vi, ai, ji); return

        idx = np.searchsorted(t, t_now, side="right") - 1
        idx = max(0, min(idx, len(t)-2))
        t0, t1 = t[idx], t[idx+1]; span = (t1 - t0)
        alpha = 0.0 if span <= 0 else (t_now - t0) / span
        xi = x[idx] + (x[idx+1]-x[idx]) * alpha
        vi = v[idx] + (v[idx+1]-v[idx]) * alpha
        ai = a[idx] + (a[idx+1]-a[idx]) * alpha
        ji = j[idx] + (j[idx+1]-j[idx]) * alpha
        self._place_marker(xi); self._update_labels(t_now, xi, vi, ai, ji)

    def _place_marker(self, x_mm):
        x0 = 40; x1 = self.canvas_w - 40
        px = x0 + (x_mm / self._xmax_canvas) * (x1 - x0) if self._xmax_canvas > 0 else x0
        y = self.canvas_h // 2; r = self.marker_r
        self.canvas.coords(self.marker_id, px - r, y - r, px + r, y + r)

    def _update_labels(self, t, x, v, a, j):
        self.lbl_t.config(text=f"Time: {t:0.1f} s")
        self.lbl_x.config(text=f"Distance: {x:0.1f} mm")
        self.lbl_v.config(text=f"Speed: {v:0.1f} mm/s")
        self.lbl_a.config(text=f"Acceleration: {a:0.1f} mm/s²")
        self.lbl_j.config(text=f"Jerk: {j:0.1f} mm/s³")

    # -------- Export --------
    def export_csv(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")],
                                                 initialfile="perfil_scurve_dist.csv", title="Guardar CSV")
        if not file_path: return
        d = self.data
        try:
            with open(file_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f); w.writerow(["t_s","v_mm_s","a_mm_s2","j_mm_s3","x_mm"])
                for i in range(len(d["t"])):
                    w.writerow([f"{d['t'][i]:.6f}", f"{d['v'][i]:.6f}", f"{d['a'][i]:.6f}", f"{d['j'][i]:.6f}", f"{d['x'][i]:.6f}"])
            messagebox.showinfo("Éxito", f"Archivo guardado:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Error saving", str(e))

    def destroy(self):
        if self._timer is not None:
            try: self.after_cancel(self._timer)
            except Exception: pass
            self._timer = None
        super().destroy()


if __name__ == "__main__":
    try:
        app = SuiteApp()
        app.mainloop()
    except Exception as e:
        print("Error:", e)
        sys.exit(1)
