import gradio as gr
from tree import molecule_breaker

res = molecule_breaker("(N(F)(N(H(H)(H))))", 5)

def greet(name, intensity):
    return f"{res}\n"  + "Hello, " + name + "!" * int(intensity)

demo = gr.Interface(
    fn=greet,
    inputs=["text", "slider"],
    outputs=["text"],
)

demo.launch()
