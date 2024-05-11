"""A very barebones copy of Codegebra for use by the model"""
import math
import re

pattern = r'-?[0-9.]+\s*=\s*(?:[0-9.]+\s*\*\s*)?[0-9.]+\^\(x-[0-9.]+\)(?:\s*\+\s*[0-9.]+)?'

def parseEq(eq, eq_type="="):
    eq = eq.replace(" ", "")  # Remove whitespace
    if eq_type == "=":
        # Equation mode
        right, left = eq.split(
            "=")  # Split string into left and right sides of the equal sign
        # Seperate terms on each side of the equation
        left = " ".join(left.split("+"))
        left = left.replace("-", " -")
        right = " ".join(right.split("+"))
        right = right.replace("-", " -")
        right = right.split()
        left = left.split()
        return simplifyTerms(right), simplifyTerms(left)
    elif eq_type == ">":
        # Inequality mode
        right, left = eq.split(
            ">")  # Split string into left and right sides of the equal sign
        # Seperate terms on each side of the equation
        left = " ".join(left.split("+"))
        left = left.replace("-", " -")
        right = " ".join(right.split("+"))
        right = right.replace("-", " -")
        right = right.split()
        left = left.split()
        return simplifyTerms(right), simplifyTerms(left)
    elif eq_type == "<":
        # Inequality mode
        right, left = eq.split(
            "<")  # Split string into left and right sides of the equal sign
        # Seperate terms on each side of the equation
        left = " ".join(left.split("+"))
        left = left.replace("-", " -")
        right = " ".join(right.split("+"))
        right = right.replace("-", " -")
        right = right.split()
        left = left.split()
        return simplifyTerms(right), simplifyTerms(left)
    elif eq_type == "≥":
        # Inequality mode
        right, left = eq.split(
            "≥")  # Split string into left and right sides of the equal sign
        # Seperate terms on each side of the equation
        left = " ".join(left.split("+"))
        left = left.replace("-", " -")
        right = " ".join(right.split("+"))
        right = right.replace("-", " -")
        right = right.split()
        left = left.split()
        return simplifyTerms(right), simplifyTerms(left)
    elif eq_type == "≤":
        # Inequality mode
        right, left = eq.split(
            "≤")  # Split string into left and right sides of the equal sign
        # Seperate terms on each side of the equation
        left = " ".join(left.split("+"))
        left = left.replace("-", " -")
        right = " ".join(right.split("+"))
        right = right.replace("-", " -")
        right = right.split()
        left = left.split()
        return simplifyTerms(right), simplifyTerms(left)
    else:
        # Expression mode
        eq = " ".join(eq.split("+"))
        eq = eq.replace("-", " -")
        eq = eq.split()
        return eq


def quadratic_formula(a, b, c):
    discriminant = b ** 2 - (4 * a * c)
    if discriminant >= 0:
        x_1 = ((-1 * b) + math.sqrt(discriminant)) / (2 * a)
        x_2 = ((-1 * b) - math.sqrt(discriminant)) / (2 * a)
    else:
        i_coeff = math.sqrt(
            abs(discriminant))  # The coefficient of i if the solution is complex
        if i_coeff >= 0:
            # If the i coefficient is positive
            if i_coeff / (2 * a) == 1:
                # If it is one
                x_1 = f"{(-1 * b) / (2 * a)}+i"
                x_2 = f"{(-1 * b) / (2 * a)}-i"
            else:
                x_1 = f"{(-1 * b) / (2 * a)}+{i_coeff / (2 * a)}i"
                x_2 = f"{(-1 * b) / (2 * a)}-{i_coeff / (2 * a)}i"
        else:
            # If the coefficient is negative
            x_1 = f"{(-1 * b) / (2 * a)}{i_coeff / (2 * a)}i"
            x_2 = f"{(-1 * b) / (2 * a)}+{(i_coeff / (2 * a)) * -1}i"
    return x_1, x_2

def simplifyTerms(terms):
    terms = [
        item for item in terms
        if item is not None and item is not False and item is not True
    ]
    x_coeff = 0  # X coefficient(aka m in y=mx+b)
    x2_coeff = 0  # X^2 coefficient(aka A in ax^2+bx+c=0)
    # Quadratic and exponential flags
    is_quadratic = False
    is_expo = "^(x" in "".join(terms)
    const = 0  # Constant(aka b in y=mx+b)
    # Variables in the exponential equation
    y = 0
    a = 0
    b = 0
    h = 0
    k = 0
    for term in terms:
        if ("x" in term
            or "X" in term) and (term.find("x**2") == -1
                                 and term.find("x^2") == -1) and not is_expo:
            # Take out the X to get the coefficient
            term = term.replace("x", "")
            term = term.replace("X", "")
            if term == "":
                # If it doesn't have a coefficient than set it to one
                x_coeff += 1
            elif term == "-":
                # take away 1 if the coefficient is just a negative sign
                x_coeff += -1
            else:
                x_coeff += float(
                    term)  # If x has a coefficient then add it to the running tally
        elif "x**2" in term or "x^2" in term:
            is_quadratic = True
            # Take out the X^2 to get the coefficient
            term = term.replace("x**2", "")
            term = term.replace("x^2", "")
            if term == "":
                # If it doesn't have a coefficient than set it to one
                x2_coeff += 1
            elif term == "-":
                # take away 1 if the coefficient is just a negative sign
                x2_coeff += -1
            else:
                x2_coeff += float(
                    term)  # If x has a coefficient then add it to the running tally
        elif is_expo:
            if len(terms) == 1:
                y = float(term)
            elif re.findall(r'\d+\^\(', term) != []:
                term = term.replace("^(", "")
                term = term.replace("x", "")
                if len(term.split("*")) == 2:
                    a = float(term.split("*")[0])
                    b = float(term.split("*")[1])
                else:
                    b = float(term)
            elif re.findall(r'\d\)', term):
                h = float(term.replace(")", "")) * -1
            else:
                k = float(term)
            # Take out the X^2 to get the coefficient
            term = term.replace("x**2", "")
            term = term.replace("x^2", "")
        else:
            const += float(term)  # Add to the running tally for the constant
    if x_coeff == 1:
        x_coeff = ""
    else:
        x_coeff = str(x_coeff)
    # If a is 0 in a quadratic, is it really a quadratic?
    if x2_coeff == 1:
        x2_coeff = "x^2"
    elif x2_coeff == 0:
        is_quadratic = False
        x2_coeff = None
    else:
        x2_coeff = str(x2_coeff) + "x^2"
    # Determine which return values should be set to null
    if is_expo:
        return [None, None, False, None, a, b, h, k]
    if x_coeff == "0" and const != "0":
        return [None, str(const), is_quadratic, x2_coeff]
    elif x_coeff != "0" and const == "0":
        return [None, None, is_quadratic, x2_coeff]
    elif x_coeff == "0" and const == "0":
        return [None, None, is_quadratic, x2_coeff]
    else:
        return [f"{x_coeff}x", str(const), is_quadratic, x2_coeff]
def quadratic_formula(a, b, c):
    discriminant = b ** 2 - (4 * a * c)
    if discriminant >= 0:
        x_1 = ((-1 * b) + math.sqrt(discriminant)) / (2 * a)
        x_2 = ((-1 * b) - math.sqrt(discriminant)) / (2 * a)
    else:
        i_coeff = math.sqrt(
            abs(discriminant))  # The coefficient of i if the solution is complex
        if i_coeff >= 0:
            # If the i coefficient is positive
            if i_coeff / (2 * a) == 1:
                # If it is one
                x_1 = f"{(-1 * b) / (2 * a)}+i"
                x_2 = f"{(-1 * b) / (2 * a)}-i"
            else:
                x_1 = f"{(-1 * b) / (2 * a)}+{i_coeff / (2 * a)}i"
                x_2 = f"{(-1 * b) / (2 * a)}-{i_coeff / (2 * a)}i"
        else:
            # If the coefficient is negative
            x_1 = f"{(-1 * b) / (2 * a)}{i_coeff / (2 * a)}i"
            x_2 = f"{(-1 * b) / (2 * a)}+{(i_coeff / (2 * a)) * -1}i"
    return x_1, x_2


def exponential_solver(y, b, a=1, h=0, k=0):
    numerator = math.log((-1 * a) / (k - y)) - h * math.log(b) + (2 * math.pi)
    return f"{math.log((-1 * a) / (k - y)) - h * math.log(b) / math.log(b)}+{(2 * math.pi) / math.log(b)}i n"

def solve(equation):
    eq_type = "="
    if ">" in equation:
        eq_type = ">"
    elif "<" in equation:
        eq_type = "<"
    elif "≥" in equation:
        eq_type = "≥"
    elif "≤" in equation:
        eq_type = "≤"
    if "x**3" in equation or "x^3" in equation:
        return "Cannot solve polynomials with a degree above 2"
        return 0
    if ("x**2" in equation) or ("x^2" in equation):
        left, right = parseEq(equation, eq_type)
        right_cpy = right.copy()  # Copy the equation to avoid an infinite loop
        # Move terms to left side of the equal sign
        if right[0] is not None:
            if "-" in right[0]:
                # Add to both sides
                right.append(right[0][1:])
                left.append(right[0][1:])
            else:
                # Subtract from both sides
                right.append("-" + str(right[0]))
                left.append("-" + str(right[0]))
        if right[1] != None:
            if "-" in right[1]:
                # Add to both sides
                right.append(right[1][1:])
                left.append(right[1][1:])
            else:
                # Subtract from both sides
                right.append("-" + str(right[1]))
                left.append("-" + str(right[1]))
        if right[3] != None:
            if "-" in right[3]:
                # Add to both sides
                right.append(right[3][1:] + "x^2")
                left.append(right[3][1:] + "x^2")
            else:
                # Subtract from both sides
                right.append("-" + str(right[3]) + "x^2")
                left.append("-" + str(right[3]) + "x^2")
        right, left = simplifyTerms(right), simplifyTerms(left)
        # Extract A, B, and C
        left[3] = left[3].replace("x^2", "")
        if left[3] == "":
            a = 1
        else:
            a = float(left[3].replace("x^2", ""))
        if left[0] == None:
            b = 0
        else:
            b = float(left[0].replace("x", ""))
        c = float(left[1])
        x_1, x_2 = quadratic_formula(a, b, c)  # Plug into the quadratic formula
        if x_1 == x_2:
            return f"x {eq_type} {x_1}"
        else:
            return f"x {eq_type} {x_1} and x {eq_type} {x_2}"
    elif re.findall(pattern, equation) != []:
        left, right = parseEq(equation, eq_type)
        return f"x {eq_type} {exponential_solver(float(left[1]), right[5], right[4], right[6], right[7])} n ∈ ℤ(ℤ is the set of integers)"
    else:
        left, right = parseEq(equation, eq_type)
        # Goal: Get X on one side and the consts on the other side, then divide the coefficient by the constant
        if right[0] != None:
            # Isolate x on the left side of the equation
            if "-" in right[0]:
                # Add to both sides
                right.append(right[0][1:])
                left.append(right[0][1:])
            else:
                # Subtract from both sides
                right.append("-" + str(right[0]))
                left.append("-" + str(right[0]))
        if right[1] != None:
            # Isolate the constants on the right side
            if "-" in left[1]:
                # Add to both sides
                right.append(left[1][1:])
                left.append(left[1][1:])
            else:
                # Subtract from both sides
                right.append("-" + str(left[1]))
                left.append("-" + str(left[1]))
        left, right = simplifyTerms(left), simplifyTerms(right)
        if left[0] == "x":
            x_coeff = 1
        else:
            x_coeff = float(left[0].replace("x", ""))
        if eq_type != "=" and x_coeff < 0:
            if eq_type == ">":
                eq_type = "<"
            elif eq_type == "<":
                eq_type = ">"
            elif eq_type == "≥":
                eq_type = "≤"
            elif eq_type == "≤":
                eq_type = "≥"
        const = float(right[1])
        return f"x {eq_type} {const / x_coeff}"

def derivative(expression):
    global ans, ans_alt
    og_expression = expression
    expression = parseEq(expression, "EXP")
    new_expression = []
    for term in expression:
        if re.findall(r'(\d+(\*?))?x\^(\d+)', term):
            # Power rule
            term = term.split("x^")
            if term[0] == "":
                a = 1
            else:
                a = int(term[0])
            b = int(term[1])
            a = a * b
            b = b - 1
            if b == 1:
                b = "x"
            elif b == 0:
                b = ""
            else:
                b = f"x^{b}"
            if a == 1:
                a = ""
            term = f"{a}{b}"
            new_expression.append(term)
            # Special cases
        elif term == "log x" or term == "log(x)":
            new_expression.append("1/x")
        elif term == "sin x" or term == "sin(x)":
            new_expression.append("cos(x)")
        elif term == "tan x" or term == "tan(x)":
            new_expression.append("sec^2(x)")
        elif term == "sinh x" or term == "sinh(x)":
            new_expression.append("cosh(x)")
        elif term == "cosh x" or term == "cosh(x)":
            new_expression.append("sinh(x)")
        elif term == "tanh x" or term == "tanh(x)":
            new_expression.append("sech^2(x)")
        elif term == "-acos x" or term == "-acos(x)" or term == "acos x" or term == "acos(x)" or term == "asin(x)" or term == "asin x":
            new_expression.append("1/√(1-x^2)")
        elif term == "atan x" or term == "atan(x)":
            new_expression.append("1/(x^2-1)")
        elif term == "acosh x" or term == "acosh(x)":
            new_expression.append("1/√(x^2+1)*√(x^2-1)")
        elif term == "asinh x" or term == "asinh(x)":
            new_expression.append("1/√(x^2+1)")
        elif term == "atanh x" or term == "atanh(x)":
            new_expression.append("1/(1 - x^2)")
        elif term == "sec^2 x" or term == "sec^2(x)":
            new_expression.append("2*tan(x)*sec^2(x)")
        elif term == "sech^2 x" or term == "sech^2(x)":
            new_expression.append("-2*tanh(x)*sech^2(x)")
        elif term == "e^x" or term == "e**x":
            new_expression.append("e^x")
        elif term == "1/x":
            new_expression.append("-1/x^2")
        elif term == "cos x" or term == "cos(x)":
            new_expression.append("sin(x)")
        elif term == "-log x" or term == "-log(x)":
            new_expression.append("-1/x")
        elif term == "-sin x" or term == "-sin(x)":
            new_expression.append("-cos(x)")
        elif term == "-tan x" or term == "-tan(x)":
            new_expression.append("-sec^2(x)")
        elif term == "-sinh x" or term == "-sinh(x)":
            new_expression.append("-cosh(x)")
        elif term == "-cosh x" or term == "-cosh(x)":
            new_expression.append("-sinh(x)")
        elif term == "-tanh x" or term == "-tanh(x)":
            new_expression.append("-sech^2(x)")
        elif term == "-asin(x)" or term == "-asin x":
            new_expression.append("-1/√(1-x^2)")
        elif term == "-atan x" or term == "-atan(x)":
            new_expression.append("-1/(x^2+1)")
        elif term == "-acosh x" or term == "-acosh(x)":
            new_expression.append("-1/√(x^2+1)*√(x^2-1)")
        elif term == "-asinh x" or term == "-asinh(x)":
            new_expression.append("-1/√(x^2+1)")
        elif term == "-atanh x" or term == "-atanh(x)":
            new_expression.append("1/(x^2-1)")
        elif term == "-sec^2 x" or term == "-sec^2(x)":
            new_expression.append("-2*tan(x)*sec^2(x)")
        elif term == "-sech^2 x" or term == "-sech^2(x)":
            new_expression.append("2*tanh(x)*sech^2(x)")
        elif term == "-e^x" or term == "-e**x":
            new_expression.append("-e^x")
        elif term == "-1/x":
            new_expression.append("1/x^2")
        elif term == "-cos x" or term == "-cos(x)":
            new_expression.append("sin(x)")
        elif re.findall(r'\d+(\*?)x', term) != []:
            # Regular slope
            if "*" in term:
                new_expression.append(term[:-2])
            else:
                new_expression.append(term[:-1])
    new_expression_str = ""
    for term in new_expression:
        if new_expression_str == "":
            new_expression_str = term
        else:
            if term[0] == "-":
                new_expression_str = new_expression_str + term
            else:
                new_expression_str = new_expression_str + "+" + term
    ans = f"d/dx ({og_expression}) = " + new_expression_str
    return f"d/dx ({og_expression}) = " + new_expression_str


def integrate(expression):
    global ans, ans_alt
    og_expression = expression
    expression = parseEq(expression, "exp")
    new_expression = []
    for term in expression:
        if re.findall(r'(\d+(\*?))?x\^(\d+)', term):
            # Inverse power rule
            term = term.split("x^")
            if term[0] == "":
                a = 1
            else:
                a = int(term[0])
            b = int(term[1])
            a = a / (b + 1)
            b = b + 1
            if b == 1:
                b = "x"
            elif b == 0:
                b = ""
            else:
                b = f"x^{b}"
            if a == 1:
                a = ""
            term = f"{a}{b}"
            new_expression.append(term)
        elif re.match(r'(\-?)\d+x', term):
            term = term.replace("x", "")
            if term == "":
                new_expression.append("x^2/2")
            else:
                new_expression.append(f"{int(term) / 2}x^2")
        # Cases that can't be solved by the inverse power rule
        elif term == "log x" or term == "log(x)":
            new_expression.append("x(log(x)-1)")
        elif term == "sin x" or term == "sin(x)":
            new_expression.append("-cos(x)")
        elif term == "e^x" or term == "e**x":
            new_expression.append("e^x")
        elif term == "1/x":
            new_expression.append("log(x)")
        elif term == "cos x" or term == "cos(x)":
            new_expression.append("sin(x)")
        elif term == "tan x" or term == "tan(x)":
            new_expression.append("-log(cos((x))")
        elif term == "sinh x" or term == "sinh(x)":
            new_expression.append("cosh(x)")
        elif term == "cosh x" or term == "cosh(x)":
            new_expression.append("sinh(x)")
        elif term == "tanh x" or term == "tanh(x)":
            new_expression.append("log(cosh(x))")
        elif term == "acos x" or term == "acos(x)":
            new_expression.append("x*acos(x)-√(1-x^2)")
        elif term == "asin x" or term == "asin(x)":
            new_expression.append("√(1-x^2)+x*asin(x)")
        elif term == "atan x" or term == "atan(x)":
            new_expression.append("x*atan(x)-0.5*log(x^2+1)")
        elif term == "acosh x" or term == "acosh(x)":
            new_expression.append("x*asinh(x)-√(x^2+1)*√(x^2-1)")
        elif term == "asinh x" or term == "asinh(x)":
            new_expression.append("x*asinh(x)-√(x^2+1)")
        elif term == "atanh x" or term == "atanh(x)":
            new_expression.append("0.5*log(1-x^2)+x*atanh(x)")
        elif term == "sec^2 x" or term == "sec^2(x)":
            new_expression.append("tan(x)")
        elif term == "sech^2 x" or term == "sech^2(x)":
            new_expression.append("tanh(x)")
        elif term == "-acos x" or term == "-acos(x)":
            new_expression.append("√(1-x^2)-x*acos(x)")
        elif term == "-log x" or term == "-log(x)":
            new_expression.append("x-x*log(x)")
        elif term == "-sin x" or term == "-sin(x)":
            new_expression.append("cos(x)")
        elif term == "-tan x" or term == "-tan(x)":
            new_expression.append("log(cos(x))")
        elif term == "-sinh x" or term == "-sinh(x)":
            new_expression.append("-cosh(x)")
        elif term == "-cosh x" or term == "-cosh(x)":
            new_expression.append("-sinh(x)")
        elif term == "-tanh x" or term == "-tanh(x)":
            new_expression.append("-log(cos(x))")
        elif term == "-asin(x)" or term == "-asin x":
            new_expression.append("-√(1-x^2)-x*asin(x)")
        elif term == "-atan x" or term == "-atan(x)":
            new_expression.append("0.5*log(x^2+1)-x*atan(x)")
        elif term == "-acosh x" or term == "-acosh(x)":
            new_expression.append("√(x^2+1)*√(x^2-1)-x*acosh(x)")
        elif term == "-asinh x" or term == "-asinh(x)":
            new_expression.append("√(x^2+1)-x*asinh(x)")
        elif term == "-atanh x" or term == "-atanh(x)":
            new_expression.append("-0.5*log(1-x^2)-x*atanh(x)")
        elif term == "-sec^2 x" or term == "-sec^2(x)":
            new_expression.append("-tan(x)")
        elif term == "-sech^2 x" or term == "-sech^2(x)":
            new_expression.append("-tanh(x)")
        elif term == "-e^x" or term == "-e**x":
            new_expression.append("-e^x")
        elif term == "-1/x":
            new_expression.append("-log(x)")
        elif term == "-cos x" or term == "-cos(x)":
            new_expression.append("-sin(x)")
        elif re.match(r'(\-?)\d+', term):
            # Reverse slope
            new_expression.append(term + "x")
    # Combine all the terms into a single string
    new_expression_str = ""
    for term in new_expression:
        if new_expression_str == "":
            new_expression_str = term
        else:
            if term[0] == "-":
                new_expression_str = new_expression_str + term
            else:
                new_expression_str = new_expression_str + "+" + term
    return f"∫ ({og_expression}) dx = {new_expression_str} + C"