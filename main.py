from langchain_community.llms import LlamaCpp
from langchain_core.callbacks import CallbackManager, StreamingStdOutCallbackHandler
import json

from simple_codegebra import integrate, derivative, solve

callback_manager = CallbackManager([StreamingStdOutCallbackHandler()])

llm = LlamaCpp(
    model_path="./models/llama-3-8b-instruct.gguf",
    n_gpu_layers=-1,
    n_batch=2,
    f16_kv=True,  # MUST set to True, otherwise you will run into problem after a couple of calls
    n_ctx=1024,
    stop=["}", "<|eot_id|>"],
    callback_manager=callback_manager,
    verbose=True,  # Verbose is required to pass to the callback manager
)
def generate_prompt(query):
    template = f"""<|begin_of_text|><|start_header_id|>system<|end_header_id|>
    You(the assistant) are an expert JSON builder designed to assist with a wide range of tasks.
    
    Assistant is able to trigger actions for User by responding with JSON strings that contain the "query" and "tool_name" parameters. The tool will return a response that the assistant can then send back to the user.
    You have access to the following actions:
    - "solve" to solve for x in an equation or inequality(only standard form quadratics, slope-intercept form linear equations, and exponentials in standard form are supported).
        To use the solve tool, you should respond like so:
        {{"tool_name": "solve", "query": "2x^2+8x-8=21"}}
    - "differentiate" to differentiate a function. If the user provides a function, only provide what comes after f(x) (i.e. 2x+5 not f(x) = 2x+5).
        To use the differentiate tool, you should respond like so:
        {{"tool_name": "differentiate", "query": "5x^7-log(x)"}}
    - "integrate" to find the indefinite integral of a function. It is key that you omit the d/dx when providing the function.
        To use the integrate tool, you should respond like so:
        {{"tool_name": "integrate", "query": "2x-4x^3"}}
    <|eot_id|><|start_header_id|>user<|end_header_id|>
    
    What is the indefinite integral of 4x^2+log(x)?<|eot_id|><|start_header_id|>assistant<|end_header_id|>
    
    {{"tool_name": "integrate", "query": "4x^2+log(x)"}} {{"tool_response": "∫ (4x^2+log(x)) dx = 1.3333333333333333x^3+x(log(x)-1) + C"}}  The indefinite integral of 4x^2+log(x) is 1.3333333333333333x^3+x(log(x)-1) + C.<|eot_id|><|start_header_id|>user<|end_header_id|>
    
    Let the equation y=-2x^2+4 represent the path of a model rocket, where y is the height of the rocket in meters and x is the time after launch in seconds. At what timestep does the rocket reach 10 meters?<|eot_id|><|start_header_id|>assistant<|end_header_id|>
    
    {{"tool_name": "solve", "query": "10=-0.5x^2+6x"}} {{"tool_response": "x = 10.0 and x = 2.0"}} The rocket reaches a height of 10 meters 2 seconds after launch.<|eot_id|><|start_header_id|>user<|end_header_id|>
    
    {query}<|eot_id|><|start_header_id|>assistant<|end_header_id|>
    """
    return template

def process_response(response):
    # try to find json in the response
    try:
        # Extract json from model response by finding first and last brackets {}
        firstBracketIndex = response.index("{")
        lastBracketIndex = len(response) - response[::-1].index("}")
        jsonString = response[firstBracketIndex - 1:lastBracketIndex + 1]
        responseJson = json.loads(jsonString)
        if responseJson['tool_name'] == 'solve':
            return f'{{"tool_response": "{solve(responseJson["query"])}"}}'
        elif responseJson['tool_name'] == 'differentiate':
            return f'{{"tool_response": "{derivative(responseJson["query"])}"}}'
        elif responseJson['tool_name'] == 'integrate':
            return f'{{"tool_response": "{integrate(responseJson["query"])}"}}'
        else:
            return f'{{"error": "No tool found with name {responseJson["tool_name"]}"}}'
    except Exception as e:
        print(e)
    # No json match, just return response
    return response

query = input("Enter your query: ")
prompt = generate_prompt(query)
response = llm(prompt)
tool_output = process_response(response + "}")
prompt = prompt + response + "} " + tool_output
llm(prompt)