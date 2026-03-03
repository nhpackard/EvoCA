"""
sdl_worker.py — SDL2 display subprocess for EvoCA.

Launched by controls.run_with_controls() via subprocess.Popen.
All output goes to the terminal where Jupyter was started.

Usage (internal):
    python sdl_worker.py <pixel_shm> <ctrl_shm> <N> <px>

ctrl_shm layout (5 × int32)
    [0] quit      1 = exit
    [1] colormode 0/1/2
    [2] step_cnt
    [3] fps × 10
    [4] paused    0/1
"""
import sys
import ctypes
import traceback


def main():
    if len(sys.argv) < 5:
        print("EvoCA SDL: bad args", flush=True)
        sys.exit(1)

    pixel_shm_name = sys.argv[1]
    ctrl_shm_name  = sys.argv[2]
    N  = int(sys.argv[3])
    px = int(sys.argv[4])
    W, H = N * px, N * px

    print(f"EvoCA SDL: starting  N={N} px={px}", flush=True)

    import numpy as np
    from multiprocessing.shared_memory import SharedMemory
    import sdl2

    # Open shared memory
    try:
        pixel_shm = SharedMemory(name=pixel_shm_name)
        ctrl_shm  = SharedMemory(name=ctrl_shm_name)
    except Exception as e:
        print(f"EvoCA SDL: SharedMemory open failed: {e}", flush=True)
        sys.exit(1)

    pixels = np.ndarray((N * N,), dtype=np.int32, buffer=pixel_shm.buf)
    ctrl   = np.ndarray((5,),     dtype=np.int32, buffer=ctrl_shm.buf)

    COLOR_MODES = ["state", "env-food", "priv-food"]

    # ── SDL2 init ─────────────────────────────────────────────────
    if sdl2.SDL_Init(sdl2.SDL_INIT_VIDEO) != 0:
        print(f"EvoCA SDL: SDL_Init failed: {sdl2.SDL_GetError()}", flush=True)
        ctrl[0] = 1
        sys.exit(1)

    print("EvoCA SDL: SDL_Init OK", flush=True)

    window_p = sdl2.SDL_CreateWindow(
        b"EvoCA",
        sdl2.SDL_WINDOWPOS_CENTERED,
        sdl2.SDL_WINDOWPOS_CENTERED,
        W, H,
        sdl2.SDL_WINDOW_SHOWN,
    )
    if not window_p:
        print(f"EvoCA SDL: SDL_CreateWindow failed: {sdl2.SDL_GetError()}", flush=True)
        ctrl[0] = 1
        sdl2.SDL_Quit()
        sys.exit(1)

    sdl2.SDL_RaiseWindow(window_p)
    print("EvoCA SDL: window created and raised", flush=True)

    surface_p = sdl2.SDL_GetWindowSurface(window_p)
    if not surface_p:
        print(f"EvoCA SDL: SDL_GetWindowSurface failed: {sdl2.SDL_GetError()}", flush=True)
        ctrl[0] = 1
        sdl2.SDL_DestroyWindow(window_p)
        sdl2.SDL_Quit()
        sys.exit(1)

    sdl2.SDL_SetSurfaceBlendMode(surface_p, sdl2.SDL_BLENDMODE_NONE)

    # Pre-compute surface pixel buffer info (doesn't change while window is open)
    surf        = surface_p.contents
    pitch_i32   = surf.pitch // 4          # row stride in int32 units
    pixels_ptr  = ctypes.cast(surf.pixels, ctypes.POINTER(ctypes.c_int32))
    dst_flat    = np.ctypeslib.as_array(pixels_ptr, shape=(H * pitch_i32,))
    dst         = dst_flat.reshape(H, pitch_i32)

    print("EvoCA SDL: entering main loop", flush=True)

    event = sdl2.SDL_Event()

    while ctrl[0] == 0:

        # Process SDL events
        while sdl2.SDL_PollEvent(ctypes.byref(event)):
            if event.type == sdl2.SDL_QUIT:
                ctrl[0] = 1
            elif event.type == sdl2.SDL_KEYDOWN:
                k = event.key.keysym.sym
                if k in (sdl2.SDLK_q, sdl2.SDLK_ESCAPE):
                    ctrl[0] = 1

        if ctrl[0]:
            break

        # Render pixels from shared memory into SDL surface
        sdl2.SDL_LockSurface(surface_p)
        src = pixels.reshape(N, N)
        if px == 1:
            dst[:N, :N] = src
        else:
            dst[:H, :W] = np.repeat(np.repeat(src, px, axis=0), px, axis=1)
        sdl2.SDL_UnlockSurface(surface_p)

        sdl2.SDL_UpdateWindowSurface(window_p)

        # Window title — hide fps when paused to avoid misleading values
        step   = int(ctrl[2])
        mode   = COLOR_MODES[min(int(ctrl[1]), 2)]
        paused = bool(ctrl[4])
        if paused:
            title = f"EvoCA  PAUSED  t={step}  color={mode}"
            sdl2.SDL_Delay(16)      # ~60 Hz ceiling; don't spin while paused
        else:
            fps   = ctrl[3] / 10.0
            title = f"EvoCA  t={step}  fps={fps:.1f}  color={mode}"
        sdl2.SDL_SetWindowTitle(window_p, title.encode())

    print("EvoCA SDL: exiting cleanly", flush=True)
    sdl2.SDL_DestroyWindow(window_p)
    sdl2.SDL_Quit()
    pixel_shm.close()
    ctrl_shm.close()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
