package org.geoserver.web.demo.function;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.geotools.feature.NameImpl;
import org.geotools.filter.FunctionFactory;
import org.opengis.feature.type.Name;
import org.opengis.filter.capability.FunctionName;
import org.opengis.filter.expression.Expression;
import org.opengis.filter.expression.Function;
import org.opengis.filter.expression.Literal;

public class LayerFunctionFactory implements FunctionFactory {    
    
    public List<FunctionName> getFunctionNames() {
        List<FunctionName> functionList = new ArrayList<FunctionName>();
        functionList.add(LayerFilterFunction.NAME);        
        return Collections.unmodifiableList( functionList );
    }    
    public Function function(String name, List<Expression> args, Literal fallback) {
        return function(new NameImpl(name), args, fallback);
    }
    public Function function(Name name, List<Expression> args, Literal fallback) {
        if( LayerFilterFunction.NAME.getFunctionName().equals(name)){
            return new LayerFilterFunction( args, fallback );
        }
        return null; // we do not implement that function
    }
}
