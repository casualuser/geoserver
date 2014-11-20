package org.geoserver.function;

import java.util.ArrayList;
import java.util.List;

import org.geotools.filter.capability.FunctionNameImpl;
import org.geotools.util.Converters;
import org.opengis.filter.capability.FunctionName;
import org.opengis.filter.expression.Expression;
import org.opengis.filter.expression.ExpressionVisitor;
import org.opengis.filter.expression.Function;
import org.opengis.filter.expression.Literal;
import org.opengis.filter.expression.VolatileFunction;

public class IsInstanceOf implements VolatileFunction, Function {

    public static FunctionName NAME = new FunctionNameImpl("isInstanceOf", Boolean.class,
            FunctionNameImpl.parameter("class", Class.class));

    private List<Expression> parameters;

    private Literal fallback;

    public IsInstanceOf() {
        this.parameters = new ArrayList<Expression>();
        this.fallback = null;
    }

    protected IsInstanceOf(List<Expression> parameters, Literal fallback) {
        this.parameters = parameters;
        this.fallback = fallback;
        
        if (parameters == null) {
            throw new NullPointerException("parameter required");
        }
        if (parameters.size() != 1) {
            throw new IllegalArgumentException(
                    "isInstanceOf(class) requires one parameter only");
        }
    }

    @Override
    public Object evaluate(Object object) {
        return evaluate(object, Boolean.class);
    }

    @Override
    public <T> T evaluate(Object object, Class<T> context) {
        Expression clazzExpression = parameters.get(0);
        
        Class clazz = clazzExpression.evaluate(object, Class.class);
        
        boolean result = false;
        
        if(clazz != null){
            if(clazz == Object.class){
                result = true;
            } else {
                result = clazz.isAssignableFrom(object.getClass());
            }
        }

        return Converters.convert(result, context);
    }

    @Override
    public Object accept(ExpressionVisitor visitor, Object extraData) {
        return visitor.visit(this, extraData);
    }

    @Override
    public String getName() {
        return NAME.getName();
    }

    @Override
    public FunctionName getFunctionName() {
        // TODO Auto-generated method stub
        return NAME;
    }

    @Override
    public List<Expression> getParameters() {
        return parameters;
    }

    @Override
    public Literal getFallbackValue() {
        return fallback;
    }

}
